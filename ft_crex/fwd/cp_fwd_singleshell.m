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
%-CREx180530

Np = length(Sdb);

% Realigned MRI by Fieldtrip interactive method if not done yet 
for i = 1 : Np
    
    % Check if mri_real already compute   
    pmri = Sdb(i).mri;
    pdir = fileparts(pmri);
    pMreal = [pdir, filesep, 'Mreal.mat'];
    if ~exist([pdir, filesep, 'mri_real.mat'], 'file') ||...
            ~exist([pdir, filesep, 'Mreal.mat'], 'file')
        % Need to read mri file (the one who was stored in BrainVisa database)
        mri_bv = read_mri(pmri);
        mri_real = realign_mri(pdir, mri_bv);

        save([pdir, filesep, 'mri_real'], 'mri_real');
        % Keep transform matrix MRI space -> realigned space
        Mreal = mri_real.transform/mri_bv.transform;  %#ok
        save(pMreal, 'Mreal');      
    end
    Sdb(i).Mreal = pMreal;

end

% Segment realigned MRI if mri_segment not done yet + define the singleshell
% head model
% A priori, not need to make the reslice MRI but maybee it was because of the
% specific MRI file use for these first tests
for i = 1 : Np
    
    % Check if mri_seg exists
    pmri = Sdb(i).mri;
    subj = Sdb(i).subj;
    pdir = fileparts(pmri);
    pfwd = make_dir([Sdb(i).dir, filesep, 'fwd']);
    pshell = [pfwd, filesep, 'vol_shell.mat'];
    if ~exist(pshell, 'file')
        % Need to read mri file (the one who was stored in BrainVisa database)
        mri_real = loadvar([pdir, filesep, 'mri_real']);
        mri_seg = segment_mri(pdir, mri_real, subj);
        

        %-- Prepare the singleshell head model from segmented MRI_realigned
        cfg = [];
        cfg.method = 'singleshell';
        vol_shell = ft_prepare_headmodel(cfg, mri_seg); 
        
        
        save(pshell, 'vol_shell');
        
        %-- Draw figure with MEG channels + head conduction volume
        pmeg = Sdb(i).meg;
        draw = filepath_raw(pmeg);
        if ~isempty(draw)
            pcor = make_dir([Sdb(i).dir, filesep, 'coreg']);
            Sgrad = ft_read_sens(draw);
            cmeg_fwd_checkreg_fig(vol_shell, Sgrad, subj, pcor)
        else
            fprintf('Unable to find raw MEG dataset inside this directory :\n%s', pmeg);
        end
    else
        Sdb(i).shell = pshell;
    end
end


function mri_seg = segment_mri(pdir, mri_real, subj)

%-- Segment the realigned MRI
cfg = [];
cfg.output= {'brain'};
cfg.spmversion = 'spm8';
mri_seg = ft_volumesegment(cfg, mri_real);  

% Figure
mri_combine = mri_real;
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

% Realign MRI using fieldtrip
function mri_real = realign_mri(pdir, mri)

cfg = [];

% Look for fiducial coordinates that was already picked
        
% Case 1 : Fiducial coordinates inside mri file (cf. with BrainStorm)
if isfield(mri,'hdr') && isfield(mri.hdr,'fiducial')...
        && isfield(mri.hdr.fiducial,'mri')
    cfg.fiducial = mri.hdr.fiducial.mri;
else
    % Check for fid mat ("fid" structure with fields "nas", "lpa", "rpa" "zpoint"
    % TO ADD = read fid.txt containing fid info too
    pfid = [pdir, filesep, 'fid.mat'];
    isfid = exist(pfid, 'file');
    if isfid
        cfg.fiducial = loadvar(pfid);
    else
        cfg.method = 'interactive';
    end
end
mri_real = ft_volumerealign(cfg, mri);

function mri = read_mri(pmri)

%-- Read MRI as imported from BrainVisa database
if strcmp(pmri(end-1:end), 'gz')
    ptest = gunzip(pmri);
    pmri = ptest{1};
end

ft_tool('freesurfer', 'add');
% Home-made function use to prevent bug with the ft_hastoolbox function that seems occured
% when other toolbox hab been added in matlab path ?? (spm12 ?)
% >> Error using ft_hastoolbox (line 450) - the FREESURFER toolbox is not installed, see http://surfer.nmr.mgh.harvard.edu/fswiki
% >> Error in ft_read_mri (line 393) >> ft_hastoolbox('freesurfer', 1);
mri = ft_read_mri(pmri);
mri.coordsys = 'ras';