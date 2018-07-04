function cp_db_import(Sdir)
% Initialized db_ft directories and files required for MEG data analysis
% 
% Check if directories and files are present in db_ft/PROJ/(group)/SUBJ 
% according to the anatomical preprocessing found in BrainVisa database
% (db_bv/PROJ/(group)/SUBJ/...
%
% Import the required files from BRAINVISA, FREESURFER and MEG DATABASES into FIELDTRIP_DATABASE
%
%-CREx180628

% Check for option
Sdp = struct('group', [], 'meg_run', []);
Sdir = check_opt(Sdir, Sdp);

if ~exist(Sdir.db_bv, 'dir')
    fprintf('No new anatomical file importation as BRAINVISA_DATABASE path not valid\npath= %s\n', Sdir.db_bv)
    isbv = 0;
else
    isbv = 1;
end

if ~exist(Sdir.db_meg, 'dir')
    fprintf('No new MEG data importation as MEG_DATABASE path not valid\npath= %s\n', Sdir.db_bv)
    ismg = 0;
else
    ismg = 1;
end

% Create FIELDTRIP database if not set yet
Sdir.db_ft = make_dir(Sdir.db_ft);

drun = Sdir.meg_run;

% Default data path to find required file to copy inside db_ft database
% from subject directory

% Default path for data inside db_brainvisa from the subject directory
pbvm = fullfile('t1mri', 'default_acquisition', 'default_analysis',...
    'segmentation', 'mesh');

dpath = [];
% MRI data path expected in brainvisa database
dpath.mri = fullfile('t1mri', 'default_acquisition');

% Default data path to find remeshed file + name of the files to find
dpath.mesh_rem = fullfile(pbvm, 'surface_analysis', '*_*white_remeshed_hiphop.gii');
% If remeshed files are not found, look for the full resolution mesh in pbvm
% directory matching with *_*white.gii (ex. subject_01_Lwhite.gii)
dpath.mesh_full = fullfile(pbvm, '*_*white.gii');

% Default texture files (for full resolution meshes)
dpath.tex_full = fullfile(pbvm, 'surface_analysis', '*_*white_parcels_marsAtlas.gii');

% MarsAtlas volume
dpath.vol = fullfile(pbvm, 'surface_analysis', '*_parcellation.nii.gz');

% Define all datapaths from ft database
dp_bv = bv_datapaths(Sdir);

% Check for transform files + anat files + meg path

Ns = length(dp_bv);
for i = 1 : Ns
    Sbv = dp_bv(i);

    % Create subject's folder if doesn't exist yet
    pft_subj = fullfile(Sdir.db_ft, Sbv.proj, Sbv.group, Sbv.subj);
    pft_subj = make_dir(pft_subj);
    
    % Check for required subdirectories in db_ft/PROJ/(group)/SUBJ
    Sout = make_mdir(pft_subj, {'mri', 'surf', 'tex', 'vol', 'trans', 'meg'});

    if isbv
        if isempty(find_files(Sout.mri, [], 0))
            % Import subject's MRI
            copy_mri(Sbv, dpath.mri, Sout.mri);
        end
        % For surf : try to find remeshed first
        if isempty(find_files(Sout.surf, [], 0))
            % Look for *_*white_remeshed_hiphop.gii in 
            % db_bv / PROJ /(group)/SUBJ/...
            %  ...t1mri/default_acquisition/default_analysis/segmentation/mesh/surface_analysis
            copy_meshes(Sbv, dpath, Sout);
        end

        % MarsAtlas volume parcellation
        if isempty(find_files(Sout.vol, [], 0))
            copy_vol(Sbv, dpath.vol, Sout.vol);
        end
        
        if isempty(find_files(Sout.trans, 'trm', 0))
            copy_trans(Sbv, Sdir.db_fs, Sout.trans);
        end
    end
    
    if ismg   
        % Special case for MEG: preproc / epochs + run directories ?

        % Import raw MEG data in specific run directories as listed in drun or
        % without run directory if specific folder(s) not found or if drun empty
        % directly in db_ft / PROJ /(group)/SUBJ/meg/continuous/(run_*)/raw
        prmeg = fullfile(Sdir.db_meg, Sbv.proj, Sbv.group, Sbv.subj);
        copy_meg(prmeg, drun, Sout.meg)
    end

end

% Copy trm files required to apply transformation matrix to initial subject's
% MRI
function copy_trans(Sbv, pdb_fs, pout)
% Referential txt file -- contain fixed paths from previous processing on a unix
% computer and with {0} for subject variable
% 
pref = fullfile(ptool, 'atlas', 'referential.txt');
fid = fopen(pref, 'r');
ctrans = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);

ctrans = ctrans{1};
Nt = length(ctrans);

% Read each transform matrix
for i = 1 : Nt
    strans = ctrans{i};
    
    % Find subject dir to cut referential paths
    isi = strfind(strans, '/{0}/');

    isbv = ~isempty(strfind(strans, 'db_brainvisa'));  %#ok
    
    if isbv
        ptr = fullfile(Sbv.dir, strans(isi+5 : end));
    else
        ptr = fullfile(pdb_fs, Sbv.proj, Sbv.group, Sbv.subj, strans(isi+5 : end));
    end
    ptr = strrep(ptr, '{0}', '*');
    if ~isempty(dir(ptr))
        copyfile(ptr, pout);
    else
        warning('\nUnable to import required trm file\nfrom: %s \nto: %s\n', ptr, pout);
    end    
end

function copy_meg(pmeg, drun, pout)

% Case without run folders (only one set of raw data inside pmeg directely)
if isempty(drun) && ~isempty(filepath_raw(pmeg))   
    praw = make_dir([pout, filesep, 'continuous', filesep, 'raw']);
    copyfile([pmeg, filesep, '*'], praw);
    return
else
    % Try to import run folders
    drun = '*';
end

% Import specific run folders

% All folders if drun=='*' OR if drun is empty and no raw data are found in 
% pmeg directory (run directories expected instead)
if ~isempty(drun) && ischar(drun) && strcmp(drun, '*')
    % Get list of all directories
    dd = dir(pmeg);    
    dnam = {dd(:).name}';
    % All folders excepted '.' and '..'
    isrd = [dd(:).isdir]' & cellfun(@(x) ~strcmp(x(1),'.'), dnam);
    drun = dnam(isrd);
end

if ischar(drun)
    drun = {drun};
end

% Copy each directories matching with drun list if not already present in pout /
% continuous / RUN / raw
cdir = {{pmeg}, 0 ; drun, 1};
prun = make_pathlist(cdir);
Nd = length(prun);

for i = 1 : Nd
    prun_in = [prun{i}, filesep, '*'];
    [~, dini] = fileparts(prun{i});        
    if length(dini)<4 || ~strcmp(dini(1:4), 'run_')
        dout = ['run_', dini];
    else
        dout = dini;
    end
    prun_out = fullfile(pout, 'continuous', dout, 'raw');
    if ~exist(prun_out, 'dir')
        prun_out = make_dir(prun_out);
        copyfile(prun_in, prun_out);
    end
end

% Copy MRI file (should be the one file with nii or gz or mri extension inside
% db_brainvisa/PROJ/(group)/SUBJ/t1mri/default_acquisition)
function copy_mri(Sbv, pbvin, pout)
pmri = find_files(fullfile(Sbv.dir, pbvin));
if isempty(pmri)
    warning('\nMRI file not found in BRAINVISA_DATABASE for subject %s\n', Sbv.subj)
    warning('Please copy the required files inside \n%s\n', pout)
else
    copyfile(pmri{1}, pout)
end

% Copy meshes and associated tex files from BRAINVISA_DATABASE to FT_DATABASE
% For decimated hiphop meshes, find the associated hiphop texture file that is
% included inside ft_crex/atlas/hiphop138 directory
function copy_meshes(Sbv, dpath, Sout)

[pmsh, isfull] = find_meshes(Sbv.dir, dpath.mesh_rem, dpath.mesh_full);

if isempty(pmsh)
    warning('\n-----\nSubj: %s\n', Sbv.subj)
    warning('Fail to import *meshes files* from brainvisa to fieldtrip database')
    warning('Please copy the required files inside \n%s\n', Sout.surf)
else
    % Copyfile!
    copyfile(pmsh, Sout.surf);
end
% Copy the associated tex files (for full resolution meshes)
if isfull
    ptex = check_files(fullfile(Sbv.dir, dpath.tex_full), 'Texture');
else 
    % Determine which decimate tex files to associated
    dmsh = dir(pmsh);
    pfile = [dmsh(1).folder, filesep, dmsh(1).name];
    ptex = cp_hiphop_tex(pfile);
end
copyfile(ptex, Sout.tex); 
        
function [pmsh, isfull] = find_meshes(psubj, pmesh_rem, pmesh_full)
isfull = 0;
% Check for remeshed file first
pmsh = check_files(fullfile(psubj, pmesh_rem), 'Remeshed surface');
% Try with full meshes
if isempty(pmsh)
    pmshf = fullfile(psubj, pmesh_full);
    warning('Looking for full resolution meshes:\ndir: %s\n', pmshf);
    pmsh = check_files(pmshf, 'Full resolution surface');           
    if ~isempty(pmsh)
        isfull = 1;
    end                            
end
      
function pfile =  check_files(pfile, styp)
dfile = dir(pfile);
if isempty(dfile)
    pfile = [];
    warning('\n%s files not found with matching path:\n%s', styp, pfile);
end

% Copy MarsAtlas volume file from BRAINVISA_DATABASE to FIELDTRIP_DATABASE
function copy_vol(Sbv, pbvin, pout)
pvol = dir(fullfile(Sbv.dir, pbvin));
if isempty(pvol)
    warning('\nMRI file not found in BRAINVISA_DATABASE for subject %s\n', Sbv.subj)
    warning('Please copy the required files inside \n%s\n', pout)
else
    copyfile([pvol(1).folder, filesep, pvol(1).name], pout)
end

function Sdir = make_mdir(proot, cdir)
Nd = length(cdir);
Sdir = [];
for i = 1 : Nd
    cnam = cdir{i};
    Sdir.(cnam) = make_dir([proot, filesep, cnam]);
end


% Find all directories for MEG data processing
function dp = bv_datapaths(Sdir)

% Main database directory
db_dir = Sdir.db_bv;

% Add group level directory if Sdir.group is not empty
if ~isempty(Sdir.group)
    cpmeg = {db_dir, 0
        Sdir.proj, 0
        Sdir.group, 0
        Sdir.subj, 1};
    igrp = 3;
    isubj = 4;
else
    cpmeg = {db_dir, 0
        Sdir.proj, 0
        Sdir.subj, 1};
    isubj = 3;
    igrp = [];    
end
iproj = 2;

[alldp, subj, grp, prj] = define_datapaths(cpmeg, isubj, igrp, iproj);

dp = cell2struct([alldp'; subj'; grp'; prj'], {'dir', 'subj', 'group', 'proj'});

