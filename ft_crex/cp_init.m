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
for i = 1 : Ns
    dps = dp_meg(i);
    dps.trans = prep_trans(dps, Sdir);
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
    
    if i==1
        Sdb = dps;
    else
        Sdb(i) = dps;
    end
end




function dps = add_anat(dps)

dps.surf = find_files([dps.dir, filesep, 'surf']);
dps.tex = find_files([dps.dir, filesep, 'tex']);
dps.vol = find_files([dps.dir, filesep, 'vol']);
pmri = find_files([dps.dir, filesep, 'mri']);
dps.mri = pmri{1};

function dps = add_mat(dps)
% Add already processed file if exist

% Realignment transformation matrix
dps = add_fpath(dps, [dps.dir, filesep, 'mri', filesep, 'Mreal.mat'], 'Mreal');

% MRI real
dps = add_fpath(dps, [dps.dir, filesep, 'mri', filesep, 'mri_real.mat'], 'mri_real');

% MRI resl
dps = add_fpath(dps, [dps.dir, filesep, 'mri', filesep, 'mri_resl.mat'], 'mri_resl');

% Conduction volume
dps = add_fpath(dps, [dps.dir, filesep, 'fwd', filesep, 'vol_shell.mat'], 'shell');

% Sources location
dps = add_fpath(dps, [dps.dir, filesep, 'fwd', filesep, 'sources.mat'], 'sources');

% MarsAtlas
dps = add_fpath(dps, [dps.dir, filesep, 'atlas', filesep, 'marsatlas.mat'], 'atlas');



function dps = add_fpath(dps, fpath, fnam)

if exist(fpath, 'file')
    dps.(fnam) = fpath;
else
    dps.(fnam) = [];
end


% Prepare transform mat for the forward model
function ptr = prep_trans(dps, Sdir)

ptrans = make_dir([dps.dir, filesep, 'trans']);

% Expected transform mat
ptr = [ptrans, filesep, 'Mtrans_ref.mat'];
dmat = dir(ptr);
if isempty(dmat)

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
