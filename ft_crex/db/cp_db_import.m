function Sdir = cp_db_import(Sdir)
% Initialized db_ft directories and files required for MEG data analysis
% 
% Check if directories and files are present in db_ft/PROJ/(group)/SUBJ 
% according to the anatomical preprocessing found in BrainVisa database
% (db_bv/PROJ/(group)/SUBJ/...
%
% Import the required files from BRAINVISA, FREESURFER and MEG DATABASES into FIELDTRIP_DATABASE
%
%-CREx-181106

% Check for option
Sdp = struct('group', [],...
            'meg_run', [],...
            'db_bv', [],...
            'db_fs', [],...
            'db_meg', []);
Sdir = check_opt(Sdir, Sdp);

% Ask if new anat importation is to be done (new subject to be added in fieldtrip database)
impanat = strcmp(questdlg('New anatomical file importation from Brainvisa & Freesurfer databases?',...
    Sdir.proj,...
    'No', 'Yes', 'No'), 'Yes');

% Check database paths, ask for it if path indicates by hand in cp_main as Sdir.db_bv or
% Sdir.db_fs parameter is not valid
if impanat
    isbv = 1;
    % Ask for Brainvisa database
    if ~exist(Sdir.db_bv, 'dir')
        Sdir.db_bv = select_path('Brainvisa', Sdir.proj);
    end
    % Ask for Freesurfer database
    if ~exist(Sdir.db_fs, 'dir')
        Sdir.db_fs = select_path('Freesurfer', Sdir.proj);
    end    
else
    isbv = 0;
end

% MEG database is required to do sources analysis (the purpose of the connectivity pipeline...) 
% Raw MEG dataset taking too much disk space (~500 Mo/run), it will be not
% imported in Fieldtrip database

% Ask if new MEG importation is to be done (new subject to be added in fieldtrip database)
impmeg = strcmp(questdlg('New MEG data importation from MEG database?',...
    Sdir.proj,...
    'No', 'Yes', 'No'), 'Yes');
if impmeg
    ismeg = 1;
    % Check for DB_MEG path
    if ~exist(Sdir.db_meg, 'dir')
        Sdir.db_meg = select_path('MEG', Sdir.proj);
    end
else
    ismeg = 0;
end

% Ask if MEG processing is to be done (require valid DB_MEG path !)
% Check for DB_MEG path
if ~exist(Sdir.db_meg, 'dir')
    domeg = strcmp(questdlg({'The path to the MEG database is not valid.', ['Do you want to process MEG data ',...
        '(if so, an explorer window will appear to select the MEG database)?']},...
    Sdir.proj,...
    'No', 'Yes', 'No'), 'Yes');
    if domeg
        Sdir.db_meg = select_path('MEG', Sdir.proj);
    else
        Sdir.db_meg = [];
    end
end

% Check for pipeline database path (DB_FIELDTRIP)
pdb = fileparts(Sdir.db_ft);
if exist(pdb, 'dir') && ~exist(Sdir.db_ft, 'dir')
    % Create FIELDTRIP database if not set yet
    Sdir.db_ft = make_dir(Sdir.db_ft);
elseif ~exist(pdb, 'dir')
    Sdir.db_ft = select_ft_db;
end

% Run(s) to find in DB_MEG
drun = Sdir.meg_run;

% Default data path to find required file to copy inside db_ft database
% from subject directory

%------------- Default path for data inside db_brainvisa from the subject directory
% !! could change depending on BV / FS evolution 
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
%---------------
% Define all datapaths from ft database
dp_bv = bv_datapaths(Sdir);

% Check for transform files + anat files + meg path
Ns = length(dp_bv);

% Initialize waitbar
wb = waitbar(0, 'Data import...', 'name', 'Data import');
wb_custcol(wb, [0.6 0 0.8]);

for i = 1 : Ns
    Sbv = dp_bv(i);

    waitbar(i/Ns, wb, ['Data import: ', Sbv.group, ' ', Sbv.subj]);
    % Create subject's folder if doesn't exist yet
    pft_subj = fullfile(Sdir.db_ft, Sbv.proj, Sbv.group, Sbv.subj);
    pft_subj = make_dir(pft_subj);
    
    % Check for required subdirectories in db_ft/PROJ/(group)/SUBJ
    Sini = make_mdir(pft_subj, {'anat', 'meg'});
    Sout = make_mdir(Sini.anat, {'mri', 'surf', 'tex', 'vol', 'trans'});
    Sout.meg = Sini.meg;
     
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
     
    % Define the pointers to the raw dataset for source analysis
    % Find the raw data according to the specific run directories as listed in drun or
    % without run directory if specific folder(s) not found or if drun is empty
    % Prepare info on raw path in raw_relpath.mat file in db_ft / PROJ /(group)/SUBJ/meg/run_*/_info
    if ismeg
        prmeg = fullfile(Sdir.db_meg, Sbv.proj, Sbv.group, Sbv.subj);
        prep_meg(prmeg, drun, Sout.meg)
    end

end

close(wb);

% Dialog box to ask for database folder selection then file explorer to choose
% the folder
function dpath = select_path(sname, proj)
uiwait(msgbox(['Select the ', sname,' database folder (which contains the ', proj, ' folder)'],...
    'Data importation'));
dpath = uigetdir(pwd, ['Select the ', sname,' database directory']);

function dpath = select_ft_db
uiwait(msgbox('Select the Fieldtrip database folder (which will contain the analysis outputs)',...
    'Data importation'));
dpath = uigetdir(pwd, 'Select the Fieldtrip database directory');

% Copy trm files required to apply transformation matrix to initial subject's MRI
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
% Prepare new MEG directories in db_ft
function prep_meg(pmeg, drun, pout)

if ~exist(pmeg, 'dir')
    warning('MEG directory not found:\n%s\n', pmeg);
    return
end

%-- Case without run folders (only one set of raw data inside pmeg directly)

if isempty(drun) && ~isempty(filepath_raw(pmeg))   
    praw = make_dir([pout, filesep, 'run_1', filesep, '_info']);
    % No run folder (raw data directely in DB_MEG/proj/(group)/subject* folder) 
    raw_relpath = []; %#ok
    save([praw, filesep, 'raw_relpath'], 'raw_relpath');
    return
end

%-- Import specific run folders
if ischar(drun)
    drun = {drun};
end

% Try to import all run folders if drun not specify
if isempty(drun) || any(cellfun(@(x) strcmp(x, '*'), drun))
    drun = {'*'};
end

% All folders if drun=='*' OR if drun is empty and no raw data are found in 
% pmeg directory (run directories expected instead)
if strcmp(drun{1}, '*')
    % Get list of all directories
    drun = dirfold(pmeg);
end

% Search for run directories and save info on raw MEG data path in 
% db_fieldtrip/proj/(group)/subj*/meg/run*/_info 
cdir = {{pmeg}, 0 ; drun, 1};
prun = make_pathlist(cdir);
Nd = length(prun);
for i = 1 : Nd
    prun_in = prun{i};
    % Declare a directory as a run one only if raw data are found inside
    %%%% TO DO : add new possible raw data file format in filepath_raw to be
    %%%% more "open" with other data recording systems
    if ~isempty(filepath_raw(prun_in))
        [~, dini] = fileparts(prun_in);   
        dout = run_out(dini);
        prun_out = make_dir(fullfile(pout, dout, '_info'));
        % Relative path to raw MEG data
        raw_relpath = dini;  %#ok
        save([prun_out, filesep, 'raw_relpath'], 'raw_relpath');
    end
end

% Define name of run directory that will be created in DB_FT (as run_*
% directory) depending on the initial run directory name in DB_MEG
function dout = run_out(dini)
dini = regexprep(dini, '(?i:^run|^ses)([_-])', '');
dout = ['run_', dini];
% Very weird case ? if length(dini)==3 && (strcmp(dini, 'run') || strcmp(dini, 'ses'))

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
    pfile = [fileparts(pmsh), filesep, dmsh(1).name];
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
pvol = fullfile(Sbv.dir, pbvin);
dvol = dir(pvol);
if isempty(pvol)
    warning('\nMRI file not found in BRAINVISA_DATABASE for subject %s\n', Sbv.subj)
    warning('Please copy the required files inside \n%s\n', pout)
else   
    copyfile([fileparts(pvol), filesep, dvol(1).name], pout)
end

% Make several directories (cdir list) in proot directory
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

