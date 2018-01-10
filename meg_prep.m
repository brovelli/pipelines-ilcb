function meg_prep(nsub)
% Commands on Frioul
% frioul_batch -n "14,15,16" -c 1 -M "[[1,3,4,5,6,7,8]]" -m meg_prep
% NB: hilbert requires Signal Toolbox
% nsub  = Subject number
%
% Andrea  22/03/2013
ft_defaults

%--------------------------------------------------------------------------
% Directories
%--------------------------------------------------------------------------
dirs.rawdata  = '/riou/work/comco/brovelli.a/Data/Neurophy/MEG_TE/';
dirs.preproc  = '/riou/work/comco/brovelli.a/Results/MEG_TE/Preprocessing/';
dirs.beh_model = '/riou/work/comco/brovelli.a/Results/MEG_TE/Beh_Model/';

%--------------------------------------------------------------------------
% Markers
%--------------------------------------------------------------------------
markers{1}.dt         = [ 3 3 ];                   % Pre Post interval in msec
markers{2}.dt         = [ 3 3 ];
markers{3}.dt         = [ 3 3 ];
markers{4}.dt         = [ 1.5 5 ];

%--------------------------------------------------------------------------
% Subjects
%--------------------------------------------------------------------------
% Subject 1
n = 1;
subjects(n).dir            = 'S1';
subjects(n).ses_dir        = { '1' '2' '3' '4' '5' '6' };
subjects(n).ses_fnames     = { 'c,rfDC' 'c,rfDC' 'c,rfDC' 'c,rfDC' 'c,rfDC' 'c,rfDC' };
subjects(n).ses_fnames_out = { 'F1' 'L1' 'L2' 'L3' 'L4' 'F2' };
subjects(n).nevt           = [ 1 2 3 1 ];                    % 1 = stimulus ; 2 = action ; 3 = reward
subjects(n).evtname        = { 'S' 'A' 'R' 'SAR' };              % S = Stim; A = Action ; R = Reward
% Subject 2
n = 2;
subjects(n).dir            = 'S2';
subjects(n).ses_dir        = { '1' '2' '3' '4' '5' '6' };
subjects(n).ses_fnames     = { 'c,rfDC' 'c,rfDC' 'c,rfDC' 'c,rfDC' 'c,rfDC' 'c,rfDC' };
subjects(n).ses_fnames_out = { 'F1' 'L1' 'L2' 'L3' 'L4' 'F2' };
subjects(n).nevt           = [ 1 2 3 1 ];                    % 1 = stimulus ; 2 = action ; 3 = reward
subjects(n).evtname        = { 'S' 'A' 'R' 'SAR' };              % S = Stim; A = Action ; R = Reward
% Subject 3
n = 3;
subjects(n).dir            = 'S3';
subjects(n).ses_dir        = { '1' '7' '4' '5' '6' '8' };
subjects(n).ses_fnames     = { 'c,rfDC' 'c,rfDC' 'c,rfDC' 'c,rfDC' 'c,rfDC' 'c,rfDC' };
subjects(n).ses_fnames_out = { 'F1' 'L1' 'L2' 'L3' 'L4' 'F2' };
subjects(n).nevt           = [ 1 2 3 1 ];                      % 1 = stimulus ; 2 = action ; 3 = reward
subjects(n).evtname        = { 'S' 'A' 'R' 'SAR' };              % S = Stim; A = Action ; R = Reward
% Subject 4
n = 4;
subjects(n).dir            = 'S4';
subjects(n).ses_dir        = { '1' '2' '3' '4' '5' '6' };
subjects(n).ses_fnames     = { 'c,rfDC' 'c,rfDC' 'c,rfDC' 'c,rfDC' 'c,rfDC' 'c,rfDC' };
subjects(n).ses_fnames_out = { 'F1' 'L1' 'L2' 'L3' 'L4' 'F2' };
subjects(n).nevt           = [ 1 2 3 1 ];                    % 1 = stimulus ; 2 = action ; 3 = reward
subjects(n).evtname        = { 'S' 'A' 'R' 'SAR' };              % S = Stim; A = Action ; R = Reward
% Subject 5
n = 5;
subjects(n).dir            = 'S5';
subjects(n).ses_dir        = { '1' '2' '3' '4' '5' '6' };
subjects(n).ses_fnames     = { 'c,rfDC' 'c,rfDC' 'c,rfDC' 'c,rfDC' 'c,rfDC' 'c,rfDC' };
subjects(n).ses_fnames_out = { 'F1' 'L1' 'L2' 'L3' 'L4' 'F2' };
subjects(n).nevt           = [ 1 2 3 1 ];                    % 1 = stimulus ; 2 = action ; 3 = reward
subjects(n).evtname        = { 'S' 'A' 'R' 'SAR' };              % S = Stim; A = Action ; R = Reward
% Subject 6
n = 6;
subjects(n).dir            = 'S6';
subjects(n).ses_dir        = { '1' '2' '3' '4' '5' '6' };
subjects(n).ses_fnames     = { 'c,rfDC' 'c,rfDC' 'c,rfDC' 'c,rfDC' 'c,rfDC' 'c,rfDC' };
subjects(n).ses_fnames_out = { 'F1' 'L1' 'L2' 'L3' 'L4' 'F2' };
subjects(n).nevt           = [ 1 2 3 1 ];                    % 1 = stimulus ; 2 = action ; 3 = reward
subjects(n).evtname        = { 'S' 'A' 'R' 'SAR' };              % S = Stim; A = Action ; R = Reward
% Subject 7
n = 7;
subjects(n).dir            = 'S7';
subjects(n).ses_dir        = { '1' '2' '3' '4' '5' '6' };
subjects(n).ses_fnames     = { 'c,rfDC' 'c,rfDC' 'c,rfDC' 'c,rfDC' 'c,rfDC' 'c,rfDC' };
subjects(n).ses_fnames_out = { 'F1' 'L1' 'L2' 'L3' 'L4' 'F2' };
subjects(n).nevt           = [ 1 2 3 1 ];                    % 1 = stimulus ; 2 = action ; 3 = reward
subjects(n).evtname        = { 'S' 'A' 'R' 'SAR' };              % S = Stim; A = Action ; R = Reward
% Subject 8
n = 8;
subjects(n).dir            = 'S8';
subjects(n).ses_dir        = { '1' '2' '3' '4' '5' '6' };
subjects(n).ses_fnames     = { 'c,rfDC' 'c,rfDC' 'c,rfDC' 'c,rfDC' 'c,rfDC' 'c,rfDC' };
subjects(n).ses_fnames_out = { 'F1' 'L1' 'L2' 'L3' 'L4' 'F2' };
subjects(n).nevt           = [ 1 2 3 1 ];                    % 1 = stimulus ; 2 = action ; 3 = reward
subjects(n).evtname        = { 'S' 'A' 'R' 'SAR' };              % S = Stim; A = Action ; R = Reward
% Subject 9
n = 9;
subjects(n).dir            = 'S9';
subjects(n).ses_dir        = { '1' '2' '3' '4' '5' '6' };
subjects(n).ses_fnames     = { 'c,rfDC' 'c,rfDC' 'c,rfDC' 'c,rfDC' 'c,rfDC' 'c,rfDC' };
subjects(n).ses_fnames_out = { 'F1' 'L1' 'L2' 'L3' 'L4' 'F2' };
subjects(n).nevt           = [ 1 2 3 1 ];                    % 1 = stimulus ; 2 = action ; 3 = reward
subjects(n).evtname        = { 'S' 'A' 'R' 'SAR' };              % S = Stim; A = Action ; R = Reward
% Subject 10
n = 10;
subjects(n).dir            = 'S10';
subjects(n).ses_dir        = { '1' '2' '3' '4' '5' '6' };
subjects(n).ses_fnames     = { 'c,rfDC' 'c,rfDC' 'c,rfDC' 'c,rfDC' 'c,rfDC' 'c,rfDC' };
subjects(n).ses_fnames_out = { 'F1' 'L1' 'L2' 'L3' 'L4' 'F2' };
subjects(n).nevt           = [ 1 2 3 1 ];                    % 1 = stimulus ; 2 = action ; 3 = reward
subjects(n).evtname        = { 'S' 'A' 'R' 'SAR' };              % S = Stim; A = Action ; R = Reward
% Subject 11
n = 11;
subjects(n).dir            = 'S11';
subjects(n).ses_dir        = { '1' '2' '3' '4' '5' '6' };
subjects(n).ses_fnames     = { 'c,rfDC' 'c,rfDC' 'c,rfDC' 'c,rfDC' 'c,rfDC' 'c,rfDC' };
subjects(n).ses_fnames_out = { 'F1' 'L1' 'L2' 'L3' 'L4' 'F2' };
subjects(n).nevt           = [ 1 2 3 1 ];                    % 1 = stimulus ; 2 = action ; 3 = reward
subjects(n).evtname        = { 'S' 'A' 'R' 'SAR' };              % S = Stim; A = Action ; R = Rewa
% Subject 12
n = 12;
subjects(n).dir            = 'S12';
subjects(n).ses_dir        = { '1' '2' '3' '4' '5' '6' };
subjects(n).ses_fnames     = { 'c,rfDC' 'c,rfDC' 'c,rfDC' 'c,rfDC' 'c,rfDC' 'c,rfDC' };
subjects(n).ses_fnames_out = { 'F1' 'L1' 'L2' 'L3' 'L4' 'F2' };
subjects(n).nevt           = [ 1 2 3 1 ];                    % 1 = stimulus ; 2 = action ; 3 = reward
subjects(n).evtname        = { 'S' 'A' 'R' 'SAR' };              % S = Stim; A = Action ; R = Reward
% Subject 13
n = 13;
subjects(n).dir            = 'S13';
subjects(n).ses_dir        = { '1' '2' '3' '4' '5' '6' };
subjects(n).ses_fnames     = { 'c,rfDC' 'c,rfDC' 'c,rfDC' 'c,rfDC' 'c,rfDC' 'c,rfDC' };
subjects(n).ses_fnames_out = { 'F1' 'L1' 'L2' 'L3' 'L4' 'F2' };
subjects(n).nevt           = [ 1 2 3 1 ];                    % 1 = stimulus ; 2 = action ; 3 = reward
subjects(n).evtname        = { 'S' 'A' 'R' 'SAR' };              % S = Stim; A = Action ; R = Rewa

%--------------------------------------------------------------------------
%  Parameters
%--------------------------------------------------------------------------
pars.rmch        = {'-A151' '-A125'};    % Otherwise leave empty []
pars.visrej      = 0;
pars.ica_ecg     = 0;
pars.lpfilter    = 250;
pars.planar      = 0; 
pars.downsample  = 1;
pars.der         = 'no';
pars.hilbert     = struct('do', 0, 'type', 'abs', 'bpfilter', [ 60 120 ], 'fname', '_HGA');

%--------------------------------------------------------------------------
% Loop over subjects, sessions and events
%--------------------------------------------------------------------------
for nses = 1:length(subjects(nsub).ses_dir)
    for ne = 1:length(subjects(nsub).nevt)
        
        %--------------------------------------------------------------
        % Preprocessing: define trials using codes from Presentation + 512 (photodiode)
        %--------------------------------------------------------------
        cfg = [];
        % It makes sure that beh and MEG data have the same trials trialfun_meg_TE
        cfg.trialfun      = 'trialfun_meg_TE';
        cfg.beh_fname     = [ dirs.beh_model subjects(nsub).dir '/beh.mat' ];
        cfg.nses          = nses;
        cfg.nevt          = subjects(nsub).nevt(ne);
        cfg.dataset       = [ dirs.rawdata subjects(nsub).dir '/' subjects(nsub).ses_dir{nses} '/' subjects(nsub).ses_fnames{nses} ];
        cfg.trialdef.pre  = markers{ne}.dt(1);
        cfg.trialdef.post = markers{ne}.dt(2);
        [ cfg ]           = ft_definetrial(cfg);
        
        %--------------------------------------------------------------
        % Preprocessing: epoch data and band pass filter
        %--------------------------------------------------------------
        if isempty(pars.rmch)
            cfg.channel = {'MEG'};
        else
            cfg.channel = {'MEG' pars.rmch{:}};
        end
        cfg.continous  = 'yes';
        cfg.detrend    = 'yes';
        cfg.lpfilter   = 'yes';
        cfg.lpfreq     = pars.lpfilter;
        % Prewhitening the data to remove 1/f trend (careful, it
        % introduces a bias in the noise)
        cfg.derivative = pars.der;
        data           = ft_preprocessing(cfg);
        
        %--------------------------------------------------------------
        % Preprocessing: downsample
        %--------------------------------------------------------------
        if pars.downsample
            cfg            = [];
            cfg.resamplefs = 1000;
            cfg.detrend    = 'no';
            data           = ft_resampledata(cfg, data);
        end
        
        %--------------------------------------------------------------
        % Preprocessing: manual rejection of MEG sensors
        %--------------------------------------------------------------
        if pars.visrej
            cfg          = [];
            cfg.method   = 'trial';
            cfg.alim     = 5e-12;
            cfg.megscale = 1;
            data         = ft_rejectvisual(cfg, data);
        end
        
        %--------------------------------------------------------------
        % Preprocessing: ICA to eliminate ECG
        %--------------------------------------------------------------
        if pars.ica_ecg
            % Downsample data before ICA
            data_tmp       = data;
            cfg            = [];
            cfg.resamplefs = 150;
            cfg.detrend    = 'no';
            data_tmp       = ft_resampledata(cfg, data_tmp);
            % Run ICA
            cfg            = [];
            cfg.method     = 'runica';
            % cfg.runica.pca = 50;
            comp           = ft_componentanalysis(cfg, data_tmp);
            % Plot components
            cfg = [];
            cfg.component  = [1:10];
            cfg.layout     = '4D248.lay';
            cfg.comment    = 'no';
            ft_topoplotIC(cfg, comp);
            % Plot loading for each component
            cfg            = [];
            cfg.channel    = {comp.label{1:5}};
            cfg.layout     = '4D248.lay';
            ft_databrowser(cfg, comp)
            % Decompose the original data as it was prior to downsampling to 150Hz
            cfg            = [];
            cfg.method     = 'runica';
            % cfg.runica.pca = 50;
            cfg.unmixing   = pinv(comp.topo);
            cfg.topolabel  = comp.topolabel;
            comp_orig      = ft_componentanalysis(cfg, data_ar);
            % Remove ECG component and reconstruct signals
            pause
            cfg            = [];
            cfg.component  = 1;
            data_prep      = ft_rejectcomponent(cfg, comp_orig);
        else
            data_prep = data;
        end
        
        %--------------------------------------------------------------
        % Preprocessing: Hilbert transform on band passed signal
        %--------------------------------------------------------------
        if pars.hilbert.do
            cfg = [];
            cfg.bpfilter      = 'yes';
            cfg.bpfreq        = pars.hilbert.bpfilter;
            cfg.hilbert       = pars.hilbert.type;
            data_prep_hilbert = ft_preprocessing(cfg,data_prep);
        end
        
        %--------------------------------------------------------------
        % Save preprocessed data and clear tmp files
        %--------------------------------------------------------------
        cd( [ dirs.preproc subjects(nsub).dir '/' ] )
        if pars.hilbert.do
            save( [ subjects(nsub).ses_fnames_out{nses} '_' subjects(nsub).evtname{ne} pars.hilbert.fname ], 'data_prep_hilbert')
        end
        save( [ subjects(nsub).ses_fnames_out{nses} '_' subjects(nsub).evtname{ne} ], 'data_prep')
        clear data data_prep data_planar comp comp_orig cfg
        
    end
end

%     %--------------------------------------------------------------
%     % Preprocessing: compute planar gradient
%     %--------------------------------------------------------------
%     if pars.planar
%         cfg = [];
%         cfg_neighbours          = [];
%         cfg.feedback            = 'yes';
%         cfg_neighbours.method   = 'distance';
%         cfg_neighbours.template = 'bti248_neighb.mat';
%         cfg.neighbours          = ft_prepare_neighbours(cfg_neighbours, data_prep);
%         cfg.planarmethod        = 'sincos';
%         data_planar             = ft_megplanar(cfg, data_prep);
%         cfg = [];
%         data_planar = ft_combineplanar(cfg, data_planar);
%     end
