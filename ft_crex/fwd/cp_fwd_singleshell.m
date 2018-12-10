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

% Realigned MRI by Fieldtrip interactive method if not done yet 
Sdb = cp_fwd_mri_prep(Sdb);

% Number of data to process
Np = length(Sdb);

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
pmeg = dps.meg.raw{1};
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