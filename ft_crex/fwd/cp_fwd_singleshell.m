function Sdb = cp_fwd_singleshell(Sdb)
% Prepare single shell conduction volume for MEG source reconstruction 
% Input : Sdb containing the list of subjects to process (including data path
% information - see cp_init)
% 
% Step 1 : MRI realignment according to fiducial markers 
%           If fid are not yet set (mri/fid.mat or inside mri data
%           (mri.hdr.fiducial.mri), the fieldtrip interactive method is used
%       --> mri_real.mat is saved in db_ft/PROJ/SUBJ/mri directory + the
%       transform matrix Mreal that will be used for atlas coregistration
%
% Step 2 : from realigned mri, the segmentation of the brain is done using
%       ft_volumesegment with 'spm8' (to avoid a bug that appear using spm12)
%       --> mri_seg.mat + figure of segmentation in db_ft/PROJ/SUBJ/mri
%
% Step 3 : the volume of conduction is defined from the mri_seg using
%       ft_prepare_headmodel and 'singleshell' method
%       --> vol_shell.mat saved inside db_ft/PROJ/SUBJ/fwd
%
% Processing are not applied if files have been previously processed at subject
% level.
%
% % TO ADD = read fid.txt containing fid info too  
%   Figure showing fid location on MRI
%
%-CREx-180704


% Number of data to process
Np = length(Sdb);

% Realigned MRI by Fieldtrip interactive method if not done yet 
for i = 1 : Np   
    % Check if mri_real already compute       
    if isempty(Sdb(i).anat.mri_real)
        % If not, prepare MRI: read the MRI file, realigne and reslice it, and 
        % keep the final Mreal transformation matrix 
        Sdb(i) = prepare_mri(Sdb(i));
    end    
end

% Segment realigned MRI if mri_segment not done yet + define the singleshell
% head model
for i = 1 : Np
    if isempty(Sdb(i).fwd.shell)
        Sdb(i) = prepare_singleshell(Sdb(i));
    end
end

function dps = prepare_singleshell(dps)

pmri = dps.anat.mri;
subj = dps.info.subj;
pdir = fileparts(pmri);

% Create fwd directory to store vol_shell.mat file 
pfwd = make_dir([dps.dir, filesep, 'fwd']);
    
presl = dps.anat.mri_resl;
% Need to read mri realigned file
mri_resl = loadvar(presl);
mri_seg = segment_mri(pdir, mri_resl, subj);

%-- Prepare the singleshell head model from segmented MRI_realigned
cfg = [];
cfg.method = 'singleshell';
vol_shell = ft_prepare_headmodel(cfg, mri_seg);        

pshell = [pfwd, filesep, 'vol_shell.mat'];
save(pshell, 'vol_shell');

% Keep data path
dps.fwd.shell = pshell;

%-- Draw figure with MEG channels + head conduction volume
% Place figure in a separated 'coreg' directory to check for the results
pmeg = dps.meg.continuous.raw{1};
draw = filepath_raw(pmeg);
fprintf('\nCo-registration figures...')
if ~isempty(draw)
    pcor = make_dir([dps.dir, filesep, 'coreg']);
    Sgrad = ft_read_sens(draw);
    cmeg_fwd_checkreg_fig(vol_shell, Sgrad, subj, pcor)
    fprintf('\nCheck for figures in %s\n', pcor);
else    
    fprintf('Unable to find raw MEG dataset to add sensor positions\npath=%s\n', pmeg);
end
        
% Segment MRI for singleshell conduction volume definition
function mri_seg = segment_mri(pdir, mri_resl, subj)
  
%-- Segment the realigned MRI
cfg = [];
cfg.output= {'brain'};
cfg.spmversion = 'spm8';
mri_seg = ft_volumesegment(cfg, mri_resl);  

% Figure
mri_combine = mri_resl;
mri_combine.seg  = double(mri_seg.brain);
mri_combine.mask = mri_combine.seg(mri_combine.seg>0);

cfg = [];
cfg.interactive   = 'no';
cfg.funparameter  = 'seg';
cfg.funcolormap   = 'jet';
cfg.opacitylim    = [0 1.5];
cfg.maskparameter = 'mask';
cfg.location = 'center';
ft_sourceplot(cfg, mri_combine);  
export_fig([pdir, filesep, 'mri_seg_', subj, '.png'], '-m1.5')
close

% Prepare the MRI for segmentation
function dps = prepare_mri(dps)

pmri = dps.anat.mri;
pdir = fileparts(pmri);       

% Need to read mri file (the one who was stored in BrainVisa database)
% + be sure units are in mm
mri_bv = read_mri(pmri);

% Be sure about the coordsys
% mri.coordsys = 'ras';
% ft_determine_coordsys
mri_real = realign_mri(pdir, mri_bv);

% Keep transform matrix MRI space -> realigned space
Mreal = mri_real.transform/mri_bv.transform;  %#ok


% Reslice MRI before doing the segmentation
% Problem with transform matrix still not clear
cfg = [];
cfg.resolution = 1;
cfg.dim = [256 256 256]; % cf. Freesurfer dim
mri_resl = ft_volumereslice(cfg, mri_real);  %#ok

% Save data in mri dir
% Save mri_real.mat file in MRI directory
preal = [pdir, filesep, 'mri_real.mat'];
save(preal, 'mri_real'); 

% Transformation to MEG space
pMreal = [pdir, filesep, 'Mreal.mat'];
save(pMreal, 'Mreal');    

% Save mri_resl.mat file in MRI directory
presl = [pdir, filesep, 'mri_resl.mat'];
save(presl, 'mri_resl');

% Keep paths
dps.anat.mri_real = preal;
dps.anat.mri_resl = presl;
dps.anat.Mreal = pMreal;
        
% Realign MRI using fieldtrip
function mri_real = realign_mri(pdir, mri)

cfg = [];

% Set default coordsys ('nas', 'lpa', 'rpa' and 'zpoint' are expected)
cfg.coordsys = 'ctf';

% Look for fiducial coordinates that was already picked
% - Case 1 : Fiducial coordinates inside mri file (cf. with BrainStorm)
% - Case 2 : fid.mat hold fiducial coordinates in mri directory
% Otherwise, 'interactive' method is set
pfid = [pdir, filesep, 'fid.mat'];
if ~exist(pfid, 'file') 
    % Try to find fid in mri.hdr field
    if isfield(mri,'hdr') && isfield(mri.hdr,'fiducial')...
            && isfield(mri.hdr.fiducial,'mri')
        fid = mri.hdr.fiducial.mri;

        cfg.fiducial = fid;
        
        % Save fid
        save(pfid, 'fid');
    else
        cfg.method = 'interactive';
    end
else
    % Check for fid mat ("fid" structure with fields "nas", "lpa", "rpa"
    % "zpoint")
    cfg.fiducial = loadvar(pfid);
end

mri_real = ft_volumerealign(cfg, mri);

if strcmp(mri_real.cfg.method, 'interactive')
    fid = mri_real.cfg.fiducial; %#ok
    save(pfid, 'fid');
end

% Be sure mri_real voxels are homogenous isotropic (cf. Fardin Afdideh suggestion)
% for segmentation to be done by segment_mri
% "Segmentation works properly when the voxels of the anatomical images are
% homogenous isotropic, i.e., have a uniform thickness for each slice"
% % if sum(diff(mri_real.dim)) ~= 0
%     cfg = [];
%     cfg.resolution = 1;
%     cfg.dim = [256 256 256]; % cf. Freesurfer dim
%     mri_resl = ft_volumereslice(cfg, mri_real);
% % end
% 
% Mtr = mri_resl.transform/mri_real.transform; 
