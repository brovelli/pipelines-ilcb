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
%-- Perform source reconstruction with beamformer (DICS)
%
%_________________________ INPUTS

% Distance between sources in subcortical regions
dsources_mm = 5;

%----------  Data paths

% The file structure of the database is set here as proposed by Andrea 
% See database_part directory for illustration

pdb = ['D:', filesep, 'databases_pipelines'];

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

% Run directories to find raw data in db_meg
%	- a list of directories to import 
%       ex. : - {'1', '3'} 
%             - or with jocker 'run_*' => all directory names matching with 'run_*'
%	- all directories '*'
%   - empty [] => no run directories
Sdir.meg_run = {'1', '6'}; %, '6'};

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
Sdir.subj = 'subject_*';

%----------  MEG data processing options
mopt = [];
%- Filtering continuous data
% Filter type: 'lp' (low-pass), 'hp' (high-pass), 'bp' (band-pass) 
% or empty '' for no filtering
mopt.continuous.filt.type = 'bp';
% Cut-off frequency in Hz (a unique number for 'hp' and 'lp', a vector of 2 number for 'bp'
% (ex. type = 'bp'; fc = [0.5 400];)
mopt.continuous.filt.fc = [0.5 400];

%- ICA component rejection 
% Do ICA and rejection: reject = 1; or don't: reject = 0;
mopt.continuous.ica.reject = 1;
% Number of components (a number ex. 40, or 'all')
mopt.continuous.ica.numcomp = 'all'; 

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
mopt.epoched.condition_dts = {'S', [2.000 1.500]
                                'A', [2.000 1.500]};

%%% TO ADD
% Baseline intervals for each conditions (cf. for study with priming, 
% default being the part of each epoch that is in the negative time part ( t < 0s)

%- Resampling of epoching data - empty to keep raw data
% Ex. resample_fs = 1000; (in Hz)
mopt.epoched.resample_fs = [];

%- Type of trial association regarding to the condition 
% The conditions could represent a time shift during a same
% experimental trial (ex. one trial = stimulus + action + response). In
% this case, removing an artefacted trial for one condition implies to
% remove the associated trials for the other conditions (= trial at the same
% position number) => rm_trials_cond = 'same'
% If the trials of the different conditions are independant, one bad trial
% list is defined and used to remove artefacted trials for each condition.
mopt.epoched.rm_trials_cond = 'same';
% ! 'same' assumes an identical number of trials per condition 

%----------  Options for source estmation by a beamforming method (LCMV or DICS)
sopt = [];
% Flag to indicate if data are to be concatenate across runs
sopt.concat_run = 0;
sopt.method = 'dics';
% Time-frequency analysis parameters for DICS
% cell with each row = one specific analysis
% At each row, 4 parameters are given:
% 1st column: analysis name as a string (ex. 'beta-gamma')
% 2nd column: frequency of interest (FOI) in Hz (ex. 30)
% 3rd column: width of the frequency smoothing in Hz (frequency resolution) (ex. 10)
% 4th column: with of the time smoothing in s (time resolution) (ex. 0.200)
% => {name, foi, Df, Dt}
sopt.dics.freq_param = {'hga', 90, 30, 0.200};

% Conditions parameters for DICS
% cell with one row per specific input data to analyse with specific time of interest
% 1st column: condition name of the input data  
% 2nd column: output data name
% 3rd column: times of interest in s as a vector of the central times of the slidding frequency 
%   analysis windows - must be compatible with the epoched condition_dts
%   parameters (with the same or shorter time boundaries)
% 4th column: normalisation indication
% => {cond_in, cond_out, toi, norm_bsl}
sopt.dics.cond_param =  {'S', 'S',  -0.700 : 0.005 : 1.500 % Stim
                        'A', 'A',   -1.000 : 0.005 : 1.000 % Action
                        'S', 'BL',  -0.500 : 0.005 : -0.100 };% Baseline

% Normalisation
% The normalisation parameters norm_param are defined in a cell with one row per 
% condition to be normalized (except the condition that is used as baseline), with:
% 1st column = cond_out name to be normalized 
% 2nd column = ¤ associated time interval to consider as baseline in second (vector of number 1x2) ex. [-0.5 0]
%              ¤ or, the cond_out name to be taken as baseline ex. 'BL'
% => {cond_out, bsl_interval OR bsl_name}
% ex. norm_param = {'S', [-0.5 0]
%                   'A', 'BL'};
%
% Special cases:
%-- For dependant conditions only (conditions = time shifts during the same trial)
% the same cond_out used as baseline can be specified by setting norm_param with
% the baseline cond_out name directly (one of the 2nd column of cond_param)
% (ex. norm_param = 'BL');
%
%-- For independant conditions, the baseline should correspond to a time interval
% outside the stimulation (ex. before the trigger)
% The baseline interval can be specified for each condition with a cell with one
% row per condition, and
%   1st col. = cond_out name 
%   2nd col. = time interval to consider as baseline [vector of number 1 x 2] in
%   seconds
% Ex. norm_param = {'S', [-0.5 0]
%                   'A', [-0.5 0]
%                  ...}
% To define a baseline interval from the trial beginning to 0 s for all conditions,
% norm_param can be let empty (norm_param = [])
%
%-- For mixed condition type (dependant/independant), 
% For each output data (cond_out) - except the condition used as baseline

sopt.dics.norm_param = 'BL';

% Don't forget to add fieldtrip toolbox to the Matlab paths (>= Mars 2018)

%______________________ END OF INPUT

% Ask for ft_crex toolbok if not presents in current working directory
ptool = 'ft_crex';
if ~exist('ft_crex', 'dir')
    ptool = uigetdir(pwd, 'Select the ft_crex toolbox directory');
end
addpath(genpath(ptool))

cp_pref
ft_defaults

% Import all required anatomical files in FIELDTRIP_DATABASE
% Check for database paths
Sdir = cp_db_import(Sdir);

% Get the subjects list to process depending on Sdir information
Sdb = cp_init(Sdir);

% Prepare MRI for fwd modelling (ask for FID if not found in directory or 
% in the header of the MRI file)
Sdb = cp_fwd_mri_prep(Sdb);

% Prepare MEG data for source estimations (if meg data have been already imported)
Sdb = cp_meg_prep(Sdb, mopt);

% Prepare single shell volume if not done yet
Sdb = cp_fwd_singleshell(Sdb);

% Prepare atlas files (surf + vol)
Sdb = cp_marsatlas(Sdb);

% Check for coregistration, correct the left/right inversion if necessary
Sdb = cp_coregistration(Sdb);

% Define sources (for subcortical and cortical parcels) 
% The verification of the co-registration must be done beforehand so 
% that the normals of the cortical sources are oriented well towards the outside 
% of the surface.
Sdb = cp_fwd_sources(Sdb, dsources_mm);

% Prepare source model including leadfield and head model => for the leadfield
% computation, we need to have the final CHANNEL selection depending on MEG data
% processing and of the run (cf. if not the same channel selection / run)
Sdb = cp_fwd_leadfield(Sdb);  

% Beamforming - only DICS for now
Sdb = cp_inv_beamforming(Sdb, sopt);

