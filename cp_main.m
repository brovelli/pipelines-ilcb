% MEG_PIPELINE project
% First tests before optimization
%
%-- Initialize the db_fieldtrip database: 
%   data importation from the brainvisa, freesurfer and meg database
%   Do nothing if data previously imported or if path to brainvisa database is
%   not valid
%
%-- Prepare the anatomical and atlas files:
%   ask once for fiducials by fieldtrip interactive method if no previously saved fid.mat file
%   is found or no fid found in the raw mri file (mri.hdr.fiducial.mri)
%   
%-- Prepare the forward model:
%   coregistration, singleshell conduction volume, subcortical and cortical
%   sources with associated atlas labels
% 
%-- Prepare MEG data: TO DO
%
%-- Compute the leadfield:
%   depending on the MEG channel selection from preprocessed MEG data
%
%

%_________INPUTS

% Distance between sources in subcortical regions
dsources_mm = 5;

%----------  Data paths

% The file structure of the database is set here as proposed by Andrea 
% See database_part directory for illustration
% At term, we should initialize all data paths with a function
% ex. : Sdir = cp_init(pdb, proj, subj);
% + We should have a copy of the transformation files inside db_ft/trans
% directories (no more need of db_brainvisa and db_freesurfer paths)

pdb = ['F:', filesep, 'db_pip_decim_test'];

Sdir = [];

% Paths to database directories

% Database directory names
% db_brainvisa and db_freesurfer are initially required to initialize the data
% importation into the FIELDTRIP_DATABASE
% If data importation was previously done, the BV and FS DATABASE paths will not
% be used.
% If the FIELDTRIP_DB doesn't exist yet, it will be created during importation
Sdir.db_bv = [pdb, filesep, 'db_brainvisa'];
Sdir.db_fs = [pdb, filesep, 'db_freesurfer'];
Sdir.db_ft = [pdb, filesep, 'db_fieldtrip'];
Sdir.db_meg = [pdb, filesep, 'db_meg'];

% Run directories to copy inside db_ft/PROJ/(group)/SUBJ/meg directory
%	- a list of directories to import 
%       ex. : - {'1', '3'} 
%             - or with jocker 'run_*' => all directory names matching with 'run_*'
%	- all directories '*'
%   - empty [] => no run directories
Sdir.meg_run = {'1', '6'};

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

%___________________ END OF INPUT


addpath(genpath('ft_crex'))
cp_pref
ft_defaults

% Import all required files in FIELDTRIP_DATABASE
% if not done yet
% = files from anatomical preprocessing in Brainvisa 
% + trm files
% + raw MEG data
% If Sdir.db_bv is not a valid path to reach brainvisa database,
cp_db_import(Sdir);

% Get the subjects list to process depending on Sdir information
Sdb = cp_init(Sdir);

% Prepare single shell volume if not yet done
Sdb = cp_fwd_singleshell(Sdb);

% Prepare atlas files (surf + vol)
Sdb = cp_marsatlas(Sdb);

% Add sources (for subcortical and cortical parcels) 
Sdb = cp_fwd_sources(Sdb, dsources_mm);

%%% TO DO: meg preprocessing before leadfield computation
% Sdb = cp_meg_preproc(Sdb, popt);

% Prepare source model including leadfield and head model => for the leadfield
% computation, we need to have the final CHANNEL selection depending on MEG data
% processing
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

