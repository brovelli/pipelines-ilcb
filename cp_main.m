% MEG_PIPELINE project
% First tests before optimization
%
%-- Initialize the db_fieldtrip database: 
%   data importation from the brainvisa, freesurfer and meg database
%   Do nothing if data previously imported or if path to brainvisa database is
%   not valid
%
%-- Prepare MEG data: bad channels, components from ICA, strong
%   artefacts and bad trials identification 
%
%-- Prepare the anatomical and atlas files:
%   ask once for fiducials by fieldtrip interactive method if no previously saved fid.mat file
%   is found or no fid found in the raw mri file (mri.hdr.fiducial.mri)
%   
%-- Prepare the forward model:
%   coregistration, singleshell conduction volume, subcortical and cortical
%   sources with associated atlas labels
%
%-- Compute the leadfield:
%   depending on the MEG channel selection from preprocessed MEG data
%

%_________INPUTS

% Distance between sources in subcortical regions
dsources_mm = 5;

%----------  Data paths

% The file structure of the database is set here as proposed by Andrea 
% See database_part directory for illustration

pdb = ['F:', filesep, 'db_pip_decim_test2'];

Sdir = [];

% Paths to database directories

% Database directory names
% db_brainvisa and db_freesurfer are initially required to initialize the data
% importation into the FIELDTRIP_DATABASE
% If data importation was previously done, the BV and FS DATABASE paths are not
% required anymore
% If the FIELDTRIP_DB doesn't exist yet, it will be created during data importation
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
Sdir.meg_run = {'1'};

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
Sdir.subj = 'subject_01';

%------------ MEG data processing options
mopt = [];
%- Filtering continuous data
mopt.continuous.filt.type = 'bp';
mopt.continuous.filt.fc = [0.5 250];

%- ICA component rejection 
mopt.continuous.ica.reject = 1;
mopt.continuous.ica.numcomp = 'all'; % a number or 'all'

% Remove the same bad channels selection per subject for all run ('same') 
% or make a channel selection for each run ('each')
mopt.continuous.rm_sens_run = 'same';

%- Functions to define trials / conditions -- project specific
% Trigfun make the link between markers values (/type) and condition names
% Trialfun allows to sort trials depending on response (to define correct
% trial for example). If trialfun is empty, the default ft_trialfun_general
% is used.
mopt.epoched.trigfun = 'trigfun_te';
mopt.epoched.trialfun = 'trialfun_te';

%- Epoching condition and associated epoching intervals (dt_s) cell of dimension
% number_of_condition x 2 with first column = condition names (as defined in
% trigfun) and second column = associated epoching intervals as a   
% [prestim postim] vector, in second with positive values 
mopt.epoched.condition_dts = {'S', [3 3]
                                'A', [3 3]
                                'R', [3 3]
                                'SAR', [1.5 5]};

%%% TO ADD
% Baseline intervals for each conditions (cf. for study with priming, 
% default being the part of each epoch that is in the negative time part ( t < 0s)
% If the number of condition is > 1 and bsl_dt_s is a 1x2 vector, the same
% baseline interval will be used for all conditions.
% If bsl_dt_s is set to empty, no baseline correction/ normalisation will
% be done.
% Time interval relative to the epoching time vector (ex. [-3 -2] s)
% mopt.epoched.bsl_dt_s = 

%- Resampling of epoching data
mopt.epoched.resample_fs = 1000; %%%%%% []

% Remove the same bad trials selection per subject for all conditions ('same')
% or define a bad trials selection for each condition ('each')
mopt.epoched.rm_trials_cond = 'same';

% Don't forget to add fieldtrip toolbox (Mars 2018)

%___________________ END OF INPUT
% Ask for ft_crex toolbok if not presents in current working directory
ptool = 'ft_crex';
if ~exist('ft_crex', 'dir')
    ptool = uigetdir(pwd, 'Select the ft_crex toolbox directory');
end
addpath(genpath(ptool))

cp_pref
ft_defaults

% Import all required files in FIELDTRIP_DATABASE
% if not done yet = files from anatomical preprocessing in Brainvisa/Freesurfer 
% + raw MEG data
% If Sdir.db_bv is not a valid path to reach brainvisa database, data importation
% is not done
cp_db_import(Sdir);

% Get the subjects list to process depending on Sdir information
Sdb = cp_init(Sdir);

% % Prepare MRI for fwd modelling (ask for FID if not found in directory or 
% % in the header of the MRI file)
% Sdb = cp_fwd_mri_prep(Sdb);

% Prepare MEG data for source estimations
Sdb = cp_meg_prep(Sdb, mopt);
% 
% % Prepare single shell volume if not done yet
% % Ask for fid information
% Sdb = cp_fwd_singleshell(Sdb);
% 
% % Prepare atlas files (surf + vol)
% Sdb = cp_marsatlas(Sdb);
% 
% % Define sources (for subcortical and cortical parcels) 
% Sdb = cp_fwd_sources(Sdb, dsources_mm);
% 
% % Prepare source model including leadfield and head model => for the leadfield
% % computation, we need to have the final CHANNEL selection depending on MEG data
% % processing and of the run (cf. if not the same channel selection / run)
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
% % sourceAll  = ft_sourceanalysis(cfg, cleanTrials.S);

