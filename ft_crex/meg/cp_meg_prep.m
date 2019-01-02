function Sdb = cp_meg_prep(Sdb, opt)
% M/SEEG preprocessing to prepare data for source estimation:
%  1 - preprocessing on continuous dataset:
%       * filtering (hp, lp or bp filter) by ft_preprocessing
%       * removing of bad channels
%       * removing of strong artefacts (required for ICA) 
%       * reject ICA component(s)
% 2 - epoching:
%       * removing of bad trials
%       * downsampling
%
% All steps are optional and configurable in the opt structure:
%
%-- for data processing on continuous dataset:
%
% - opt.continuous.filt: filtering option with fields:
%   .type: filter type ('hp', 'lp' or 'bp') 
%       [ default: [] (empty == no filtering) ]
%	.fc: cut-off frequency(ies) in Hz - one value for 'lp' and 'hp' type,
%       2 values for 'bp' [default: []]
%
% - opt.continuous.ica: ICA parameters with fields:
%   .reject: flag to indicate if ICA rejection is to be done (when reject==1) 
%       [ default: 0 (no ICA)]
%   .numcomp: number of components
%       [ default: 'all' ]
%
% - opt.continuous.rm_sens_run: 'same' or 'each': string that indicate if
% the same bad channels are to be considered for all the runs of the 
% subject recording session or if different selections of bad channels are
% to be set for each run [ default: 'each' - bad channel selection is asked for
%   each run ]
%
%-- epoching parameters: 
%
% - opt.epoched.trigfun: project-specific function to define markers values
%       and associated condition names for epoching
%       [ see ft_crex/_specproj/trigfun_te.m for example ]
%
% - opt.epoched.trialfun: project-specific function to define trials
% according to response triggers values 
%       [ see ft_crex/_specproj/trialfun_te.m for example - default is the
%       'ft_trialfun_general' ]
%
% - opt.epoched.conditions: to specify which conditions as defined in trigfun
%       functions are to be processed, and in the same order as defined in 
%       the conditions cell. Ex.: {'A', 'S', 'R', 'SAR'} 
%       [ default: [] (empty == all conditions in trigfun function are to 
%       process) ]
%
% - opt.epoched.dt_s: [prestim postim] durations relative to the trigger in
%   seconds. [ default: [0.5 1] s ]
%   * If several conditions were specified in opt.epoched.conditions but only
%   one dt_s interval, the same epoching interval is applied for all the
%   conditions.
%   * For epoching of different lengths according to the condition,
%   the opt.epoched.conditions must be set with the associated opt.epoched.dt_s
%   holding the epoching intervals par conditions at each row (
%   ex.: [3 3; 3 3; 3 3; 1.5 5] associated with {'A', 'S', 'R', 'SAR'}
%
% - opt.epoched.rm_trials_cond: 'same' or 'each' to indicate if bad trials
% selection that was done for the first condition is to be applied to all the
% other conditions (for the case when conditions == different epoching
% versions/durations with minor time shifting)
%   [ default: 'each' - bad trials are not the same per condition ]
%
% - opt.resample_fs: downsampling frequency in Hz 
%       [ default: [] (no remsampling)
%
% When new MEG data are found during database initialisation by cp_init.m,
% the preprocessing pipeline is executed according to opt otpions.
%
% A set of figures is generated in each dataset directory (in db_fieldtrip/
% PROJ/(group)/SUBJ/meg/contnuous/prep/(filt* or no_filt)/_preproc_fig).
% Theses figures helps for bad channels identification (see cmeg_cleanup_fig).
% At the end of creating figures for all new data sets (which can be long),
% the bad channels selection is requested in the command window.
% --> this create the file 'rm_sens.txt' in data prep directory containing the bad
% channels list as well as the 'preproc.mat' file that holds the effective data
% preprocessing options (channels, ICA components and strong artefacts
% removing). It is possible to add/modify channels to remove by
% adding/removing channels in rm_sens.txt file. According to
% opt.continuous.rm_sens_run option, the update of the channel to remove
% will be done.
% , rm_sens, rm_comp, 
%
%-- TO DO:  - GUI for artefact identification
% (to be added: - opt.epoched.tshift_s: shift time for triggers)
%
%-CREx180726
imeg = is_meg(Sdb);
if ~any(imeg)
    return
else
    Sdbm = Sdb(imeg);
end
%%% TO DO: case where no meg data found
% Default options
opt = check_opt(opt, struct('continuous', struct('filt', [], 'ica', [], 'rm_sens_run', 'each'),...
                            'epoched', []));
%-- Default for continuous dataset
%- Filtering
opt.continuous.filt = check_opt(opt.continuous.filt,...
                        struct('type', [], 'fc', []));
%- ICA
opt.continuous.ica = check_opt(opt.continuous.ica,...
                        struct('reject', 0, 'numcomp', 'all'));

%-- Default for epoching
opt.epoched = check_opt(opt.epoched, struct('trigfun', [],...
                                    'trialfun', [],...
                                    'condition_dts', [],...
                                    'rm_trials_cond', 'each',...
                                    'resample_fs', []));
                                
if isempty(opt.epoched.trigfun)
    warning('MEG preprocessing:');
    warning('opt.epoched.trigfun is required to know the marker values to consider for epoching');
    error('Abort processing')
end

if isempty(opt.epoched.condition_dts)
    warning('MEG preprocessing:');
    warning(['opt.epoched.condition_dt is required to know the ',...
        'conditions to process and the associated epoching intervals']);
    error('Abort processing')
end

%-- Fixed option for continuous data preprocessing
popt = [];
popt.type = 'hp';
popt.fc = 0.5;   
popt.res_fs = 400;
opt.preproc = popt;

%-- Check for HP filtering if ICA is required
fica = opt.continuous.ica.reject;
fopt = opt.continuous.filt;

fc_hp = 0.5;
if fica
    fopt = check_ica_filt_opt(fopt, fc_hp);
    opt.continuous.filt = fopt;
end

%-- Option for epoching (unfold condition_dts table)
opt.epoched.conditions = opt.epoched.condition_dts(:, 1);
opt.epoched.dt_s = vertcat(opt.epoched.condition_dts{:, 2});


%-- Do MEG preprocessing pipeline   
Sdbm = prep_pipeline(Sdbm, opt);

% Final review
%---- Confirm a last time all the preprocessing parameters for all subjects 
% (if any change in parameters, and ICA is required, check 
[Sdbm, ischg] = cp_meg_review_gui(Sdbm, opt);

% If any change and ICA is required, check for new preproc option
% compatibility with the one used before ICA computation
if fica && ischg
    Sdbm = prep_pipeline(Sdbm, opt);
end

%-- Do the preprocessing with the requested filtering option and the cleaning
% paramters
Sdbm = cp_meg_epoching(Sdbm, opt);

Sdb(imeg) = Sdbm;

%%%% TO DO: add the possibility to do beamforming on continuous data (cf.
%%%% resting state / etc) --> if opt.epoched is empty / or
%%%% opt.epoched.condition_dts / or an option opt.bmf.continuous = 'yes' ?

%-- Pipeline for preprocessing MEG data (selection of bad channels, bad ICs, bad
% trials
function Sdb = prep_pipeline(Sdb, opt)

%- Initialize all param_txt files
Sdb = cp_meg_prep_init(Sdb, opt);

%- Strong artefact identification (required first before spectra
% superimposition for bad channel selection) 
sopt = opt.preproc;
sopt.rm_sens_run = opt.continuous.rm_sens_run;
Sdb = cp_meg_cleanup_art(Sdb, sopt);

%- Bad channel identification
Sdb = cp_meg_cleanup_sens(Sdb, sopt);

%- ICA for artefact cleaning
sopt.Nc = opt.continuous.ica.numcomp;
Sdb = cp_meg_cleanup_ica(Sdb, sopt);

%- Bad trials identification with semi-automatic method -- on the HP and
% resampling data version, cleaning from artefacts
Sdb = cp_meg_cleanup_trials(Sdb, opt);
