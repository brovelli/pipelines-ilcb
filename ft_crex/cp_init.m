function Sdb = cp_init(Sdir)
% Initialized data paths required for data processing
% according to the specific directories structure (see the database_part
% directory for illustration)
% Sdb will handle the list of subjects directories for processing - including
% paths to the anatomical/atlas files that are searched inside
% db_ft/PROJ/SUBJ/mri ; surf ; tex and vol subdirectories
% 
% Prepare Mtrans_ref transformation matrix according to referencial.txt file
% for atlas surf and vol realignment
%
%-CREx180530

% Check for required files in databases
Sdp = struct('iproc', [], 'group', []);
Sdir = check_opt(Sdir, Sdp);

% Define all datapaths from ft database
dp_meg = meg_datapaths(Sdir);

% Check for transform files + anat files + meg path
Ns = length(dp_meg);

% Initialize waitbar
wb = waitbar(0, 'Check for fieldtrip database...', 'name', 'Initializing data paths');
wb_custcol(wb, [0.6 0 0.8]);
for i = 1 : Ns
    dpsi = dp_meg(i);
    
    dps = [];
    dps.dir = dpsi.dir;
    dps.info = rmfield(dpsi, 'dir');
    dps.sinfo = fullfile(dpsi.proj, dpsi.group, dpsi.subj);
    
    waitbar(i/Ns, wb, ['Data paths: ', dps.sinfo]);
    
    dps.anat.trans = prep_trans(dpsi, Sdir);
    % Add mri, surf, vol and tex files paths
    % File expected with extension gii or nii or gz
    
    dps = add_anat(dps);
    
    % Add already process file path (MAT) if exists
    dps = add_mat(dps);
    
    % Expected MEG data directory
    % Check for raw data ==> create grad.mat, hs.mat and event.mat for coreg
    % figures + the raw meg datapath list + the continuous + epoched that is
    % concat
    dps.meg = cp_db_megpaths(make_dir([dps.dir, filesep, 'meg']));
    
    % Initialize fwd model mat paths
    dps.fwd.model_run = cell(dps.meg.Nrun, 1);
    
    if i==1
        Sdb = dps;
    else
        Sdb(i) = dps;
    end
end

close(wb);

% Path to anatomical files in db_fieldtrip
function dps = add_anat(dps)
dpa = [dps.dir, filesep, 'anat', filesep];
dps.anat.surf = find_files([dpa, 'surf']);
dps.anat.tex = find_files([dpa, 'tex']);
dps.anat.vol = find_files([dpa, 'vol']);
pmri = find_files([dpa, 'mri']);
dps.anat.mri = pmri{1};

function dps = add_mat(dps)
% Add already processed file if exist

%---- Anat directories
pan = [dps.dir, filesep, 'anat', filesep];

% Realignment transformation matrix
dps.anat = add_fpath(dps.anat, [pan, 'mri', filesep, 'Mreal.mat'], 'Mreal');

% MRI real
dps.anat = add_fpath(dps.anat, [pan, 'mri', filesep, 'mri_real.mat'], 'mri_real');

% MRI resl
dps.anat = add_fpath(dps.anat, [pan, 'mri', filesep, 'mri_resl.mat'], 'mri_resl');

% MarsAtlas
dps.anat = add_fpath(dps.anat, [pan, 'atlas', filesep, 'marsatlas.mat'], 'atlas');

%---- Fwd directories
dps.fwd = [];
% Conduction volume
dps.fwd = add_fpath(dps.fwd, [dps.dir, filesep, 'fwd', filesep, 'vol_shell.mat'], 'shell');

% Sources location
dps.fwd = add_fpath(dps.fwd, [dps.dir, filesep, 'fwd', filesep, 'sources.mat'], 'sources');

function dps = add_fpath(dps, fpath, fnam)

if exist(fpath, 'file')
    dps.(fnam) = fpath;
else
    dps.(fnam) = [];
end


% Prepare transform mat for the forward model
function ptr = prep_trans(dps, Sdir)

ptrans = make_dir([dps.dir, filesep, 'anat', filesep, 'trans']);

% Expected transform mat
ptr = [ptrans, filesep, 'Mtrans_ref.mat'];
if ~exist(ptrans, 'file')

    % trm files are expected in standard BV and FS databases according to
    % referential.txt file // IF TRM files are always the same from BV/FS pipeline, 
    % we should integrate the referential file inside ft_crex toolbox
    opt = [];
    opt.db_dir.bv = fullfile(Sdir.db_bv, dps.proj, dps.group, dps.subj);
    opt.db_dir.fs = fullfile(Sdir.db_fs, dps.proj, dps.group, dps.subj);
    opt.trans_dir = ptrans;
    Mtrans_ref = read_trans(opt);
    if ~isempty(Mtrans_ref)
        save(ptr, 'Mtrans_ref');
    else
        ptr = [];
    end
end

% Find all directories for MEG data processing
function dp = meg_datapaths(Sdir)

% Main database directory
db_dir = Sdir.db_ft;

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

iproc = Sdir.iproc;

if ~isempty(iproc) && max(iproc) <= length(alldp)
    dp = dp(iproc);
end
