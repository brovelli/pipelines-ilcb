% MEG_PIPELINE project
% First tests before optimization
% Coregistration, forward model and subcortical sources
% Matlab 2017b

% Distance between sources in subcortical regions
dsources_mm = 5;

%----------  Data paths

% The file structure of the database is set here as proposed by Andrea 
% See database_part directory for illustration
% At term, we should initialize all data paths with a function
% ex. : Sdir = cp_init(pdb, proj, subj);
% + We should have a copy of the transformation files inside db_ft/trans
% directories (no more need of db_brainvisa and db_freesurfer paths)

pdb = [pwd, filesep, 'databases_part_decim'];

Sdir = [];

% Paths to database directories

% Database directory names
% db_bv and db_freesurfer are initially required to find transformation matrices as
% indicated in the referential.txt file (in db_ft/PROJ/referential directory)
% FOR FURTHER VERSION: If transform files are already in db_ft/PROJ/(group)/SUBJ/trans directory, these
% databases would not be required anymore + see if transform are always the same
% from the BV/FS pipeline to provide a standard referential.txt file directly in ft_crex toolbox
Sdir.db_bv = [pdb, filesep, 'db_brainvisa'];
Sdir.db_fs = [pdb, filesep, 'db_freesurfer'];
Sdir.db_ft = [pdb, filesep, 'db_ft'];

% Project directory
Sdir.proj = 'meg_te';

% Group directories which could contain subject directories (ex. {'dys', 'nr'})
% from db_ft/PROJ directory
% Empty: no group directories - directly subjects
Sdir.group = [];

% Subject directory - could be:
% -- a unique string ('subject_01' )
% -- a string with a jocker '*' (ex. 'subj*')
%   --> all subject directories with matching names are searched inside project directory
% -- a cell of string: {'subject_01', 'subject_04'} or {'control*', 'patient*'}
Sdir.subj = 'subj*';

% Indices of subjects to process if several directories are found according to
% dir.subjects indication
% Empty [] to process all subjects matching with dir.subjects list
Sdir.iproc = [];

% Don't forget to add fieldtrip toolbox (Mars 2018)

%___________________

addpath(genpath('ft_crex'))
cp_pref
ft_defaults

% Get the subjects list to process depending on Sdir information
Sdb = cp_init(Sdir);

% Prepare single shell volume if not yet done
Sdb = cp_fwd_singleshell(Sdb);

% Prepare atlas files
Sdb = cp_marsatlas(Sdb);

% Add sources (for subcortical parcels for now) 
Sdb = cp_sources_subcort(Sdb, dsources_mm);

% Prepare source model including leadfield and head model
Sdb = cp_fwd_leadfield(Sdb);

%... To be continued...
% Next step => inverse problem
%
% Ex. from the model.mat of one subject
% % disp('LCMV source analysis')
% % cfg        =  [];
% % cfg.method = 'lcmv';
% % cfg.grid   = model.subcortical.grid;
% % cfg.headmodel    = model.subcortical.headmodel;
% % cfg.lcmv.lambda     = '5%';
% % cfg.lcmv.keepfilter = 'yes';
% % cfg.lcmv.fixedori   = 'yes'; % L'orientation est calculee via une ACP // orientation 3 dipoles
% % sourceAll  = ft_sourceanalysis(cfg, avgTrialsAll);

