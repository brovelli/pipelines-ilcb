function Sdb = cp_fwd_mri_prep(Sdb)
% Prepare MRI for forward modelling
% Input : Sdb containing the list of subjects to process (including data path
% information - see cp_init)
% 
% MRI realignment according to fiducial markers 
% If fid are not yet set (mri/fid.mat or inside mri data
% (mri.hdr.fiducial.mri), the fieldtrip interactive method is used
%  --> mri_real.mat is saved in db_ft/PROJ/SUBJ/mri directory + the
%  transform matrix Mreal that will be used for atlas coregistration
% % TO ADD = read fid.txt containing fid info too  
%   Figure showing fid location on MRI
%
%-CREx-180820


% Number of data to process
Np = length(Sdb);

% Initialize waitbar
wb = waitbar(0, 'Prepare MRI...', 'name', 'MRI preprocessing');
wb_custcol(wb, [0 0.6 0.8]);

% Realigned MRI by Fieldtrip interactive method if not done yet 
for i = 1 : Np   
	waitbar(i/Np, wb, ['MRI preparation: ', Sdb(i).sinfo]);  
    % Check if mri_real already compute       
    if isempty(Sdb(i).anat.mri_real)
        % If not, prepare MRI: read the MRI file, realigne and reslice it, and 
        % keep the final Mreal transformation matrix 
        Sdb(i) = prepare_mri(Sdb(i));
    end  
end

close(wb);
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
mri_real = realign_mri(pdir, mri_bv, dps.sinfo);

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
function mri_real = realign_mri(pdir, mri, idnam)

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
        % Inform user to proceed fiducial selection
        uiwait(msgbox({'\fontsize{12}Please select the fiducials for subject: ';...
            ['\fontsize{13}\bf ', prep_tex(idnam)]}, 'MRI realign', 'help',...
            struct('WindowStyle', 'non-modal', 'Interpreter', 'tex')));
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
