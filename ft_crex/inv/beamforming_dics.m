function beamforming_dics(bopt)
% DICS beamforming
%-- bopt.param: parameters for DICS beamforming with fields:
%   
%      .freq_param: frequency parameters for frequency analysis
%           {'name', foi_Hz, Df_Hz, Dt_s ;...}  (one row per analysis)
%      .cond_param: input/output condition name and associated times of interest 
%           {'cond_in', 'cond_out', t_beg : inc : t_end ; ...}
%      .norm_cond_out: name of the output condition to be taken as baseline for
%       normalization
%      
%-- bopt.info: info (subject name / run) to be added on figure title
%
%-- bopt.acond: available condition to check the consistency with param.cond_param
%
%-- bopt.fwd: path to the forward model containing singleshell, grid and leadfield 
%       (as output by cp_fwd_leadfield)
%
%-- bopt.data: path to the cleaned data
%
%-- bopt.dir: path of results folder
%
%-CREx-181115
 %%%% TO DO compare logz and z (why log z ?)
 %%%% TO DO check for ntap K etc...
 %%%% TO DO: external function to find the directory of analysis results
 %%%% depending on paremeters for further analysis (connectivity)
 
 %___ Default options
 dopt = struct('param', struct('freq_param', [],...
                                'cond_param', [],...
                                'norm_cond_out', []),...
                'info', '',...
                'acond', [],...
                'fwd', [],...
                'data', [],...
                'new_clean', 1,...
                'dir', pwd);

%___ Check for option
bopt = check_opt(bopt, dopt);
%%% TO DO: define standard frequency band analysis if param.freq_param is empty
% Abort if any required field is empty
if ~check_bopt(bopt)
    warning('Inclomplete parameters structure for DICS analysis... Aborting calculation...');
    return;
end

%___ Prepare dics parameters and data paths
% Check for parameters, check if computation was not done before
Sdp = prep_dics(bopt);

if isempty(Sdp)
    return;
end

% Load clean data and forward model 
%- clean dataset
cleanTrials = loadvar(bopt.data);

%- forward model (with .cortical and .subcortical)
fwd_model = loadvar(bopt.fwd);

%___ DICS beamforming for each frequency bands and condition to output
[fqana, Nf] = get_names(Sdp);
% Loop over frequency analysis
for i = 1 : Nf
    % Analysis name
    fnam = fqana{i};
    
    Sd = Sdp.(fnam);   
    [cond, Nc] = get_names(Sd);

    % Loop over conditions
    for j = 1 : Nc
        % Condition name (out)
        cdout = cond{j};
        
        % Parameters
        Sdc = Sd.(cdout);     
        
        % Trials to use
        cdin = Sdc.cond_in;
        trials = cleanTrials.(cdin);
        
        %---- Spectral analysis
        cfg = Sdc.cfg;
        cfg.channel = trials.label;
        tf_trials = ft_freqanalysis(cfg, trials);
        %%% TO DO: do taper analysis outside FT to gain time (cf. Dee scripts ?)
        
        % Keep trialinfo
        if isfield(trials, 'trialinfo')
            tf_trials.trialinfo = trials.trialinfo;
        end
        
        % Loop over time samples and dipoles        
        Spow = bmf_dics(tf_trials, fwd_model, Sdc.info);      

        % Baseline normalization
        if Sdc.norm_bsl.do
            if Sdc.isbsl
                Sbsl = bsl_calc_save(Spow, Sdc);               
                continue
            end
            if ~isempty(Sdc.norm_bsl.precomp)
                Sbsl = Sdc.norm_bsl.precomp;
            end
            Spow = bsl_norm(Spow, Sbsl);
        end 
        % Save Spow and do figures if normalized power
        save_pow(Spow, Sdc);         
    end
end
    
% Fast Andrea's Beamforming   
function Spow = bmf_dics(trials_tf, fwd_model, info)

ntr = size(trials_tf.cumtapcnt, 1);
ntap = trials_tf.cumtapcnt(1);
nchan = length(trials_tf.label);
ntime = length(trials_tf.time);

% Concatenate cortical and subcortical (to avoid the loops over time and dipoles to be
% done for each case, we proceed the two in the same loops)
% Extract inside leadfield and concatenate cortical and subcortical
lf_dip = concat_ldf(fwd_model);

% Number of dipoles
ndip = length(lf_dip);

pow_full = NaN(ntr, ndip, ntime);

for j = 1 : ntime
    % Fourier spectrum at time t (1 freq here)
    Fsp_ini = squeeze(trials_tf.fourierspctrm(:, :, 1, j)); % [(Ntapers*Ntrials) x Nchan]

    %--- Compute filter from avg trials
    % Cfa is computed using fixcsd in ft_checkdata with method 'fullfast' --> l.976-998

    %-- Cross-spectrum matrix (mean across trials and tapers)
    Fspa = Fsp_ini;
    Fspa(~isfinite(Fspa)) = 0; % ? really usefull ?)
    Fspa = transpose(Fspa);
    n   = sum(Fspa~=0,2);  % = Ntrials*Ntapers = 583
    Cfa = Fspa*Fspa'./n(1);

    % Regularization parameter
    ratio = 0.05;
    lambda = ratio * trace(Cfa)/size(Cfa,1);
    %-- Inverse cross-spectrum matrix by keeping real part only to compute real filter
    % strcmp(realfilter, 'yes')
    % the filter is computed using only the leadfield and the inverse covariance or CSD matrix
    % therefore using the real-valued part of the CSD matrix here ensures a real-valued filter
    invCfa = pinv(real(Cfa) + lambda * eye(size(Cfa)));
    %%%% TO DO: compute a filter for all conditions ? 
    %-- Compute filter at each dipole location
    filt_dip = cell(ndip, 1);
    for k = 1 : ndip
        lf = lf_dip{k};
        % compute filter
        filt = pinv(lf' * invCfa * lf) * lf' * invCfa;   % Gross eqn. 3, use PINV/SVD to cover rank deficient leadfield
        % Find the optimal orientation based on SVD if lf is in 3 directions
        % [Nchan x 3] 
        % (Fieldtrip code)
        if size(lf, 2)==3
            u = svd(real(filt) * Cfa * ctranspose(filt));
            maxpowori = u(:, 1);
            % Compute the leadfield for that orientation
            lf = lf * maxpowori;
            % Recompute the filter to only use that orientation
            filt = pinv(lf' * invCfa * lf) * lf' * invCfa;
        end
        filt_dip{k} = filt;
    end

    %--- Apply filter to each trials and compute power

    %-- Cross-spectra matrix for trials (mean across tapers)
    % Cf is computed using fixcsd in ft_checkdata with method 'full' --> l.946-961
    % Change fourier spectrum dim 
    % [(Ntrials*Ntrapers) x Nchan ]--> [Ntrials x Ntapers x Nchan]
    % [583 x 243 ] --> [53 x 11 x 243 ]
    Fsp = permute(reshape(Fsp_ini, ntap, ntr, nchan), [2 1 3]); 
    %--> [Ntrials x Ntapers x Nchan]

    %----------------------
    % ANDREA
    %----------------------
    A = reshape([filt_dip{:}], nchan, ndip)'; 
    S = zeros(ntr, ntap, ndip);  
    for itap = 1 : ntap
        S(:, itap, :) = transpose( A * transpose(squeeze(Fsp(:, itap, :))) );
    end
    % Compute power for each taper
    S = S.*conj(S);

    % Average over tapers 
    pow_full(:, :, j) = squeeze(mean(S, 2));   
end
% Re-split cortical and subcortical
Spow = split_pow(fwd_model, pow_full, trials_tf, info);

% Re-split cortical and subcortical part of the full power matrix
% Add info + mesh for further processing (dynmesh figures, connectivity...)
function Spow = split_pow(fwd_model, pow_full, trials, info)
[ntr, ~, nti] = size(pow_full);
% Initialize full pow matrix
%- cortical
insco = fwd_model.cortical.grid.inside;
ndip_co = length(insco);
powco = NaN(ntr, ndip_co, nti); 
%- subcortical
inssu = fwd_model.subcortical.grid.inside;
ndip_su = length(inssu);
powsu = NaN(ntr, ndip_su, nti); 

% Fill pow with values
iss = sum(insco) + 1;
powco(:, insco, :) = pow_full(:, 1:iss-1, :);
powsu(:, inssu, :) = pow_full(:, iss:end, :);

pow = [];
pow.cortical = powco;
pow.subcortical = powsu;

[type, Nt] = get_names(fwd_model);
time = trials.time;

Spow = [];
for i = 1 : Nt
    ctyp = type{i};
    
    fwd = fwd_model.(ctyp);
    
    % Keep mesh for source plot
    sop = [];
    sop.vertices = fwd.grid.pos;
    % tri and mom field: only for cortical sources
    if isfield(fwd.grid, 'tri')
        sop.faces = fwd.grid.tri;
        sop.norm = fwd.grid.mom';
    end
    
    Sp = [];
    Sp.info = info;
    Sp.time = time;
    Sp.pow = pow.(ctyp);
    Sp.inside = fwd.grid.inside;
    Sp.mesh = sop;
    Sp.atlas = fwd.atlas;
    % Keep trialinfo if presents
    if isfield(trials, 'trialinfo')
        Sp.trialinfo = trials.trialinfo;
    end
    Spow.(ctyp) = Sp;
end

% Concate cortical and subcortical leadfield to compute power in the same time
% loop
function ldf = concat_ldf(fwd_model)
insco = fwd_model.cortical.grid.inside;
inssu = fwd_model.subcortical.grid.inside;
ldf = [fwd_model.cortical.grid.leadfield(insco) fwd_model.subcortical.grid.leadfield(inssu)];

% Calculate the baseline correction metrics (mean and std)
% Save results for further source analysis with precomputed baseline
function Sbsl = bsl_calc_save(Spow, Sdc)

[typ, Nty] = get_names(Spow);

%- compute metrics for baseline correction
Sbsl = [];
for i = 1 : Nty
    ctyp = typ{i};
    pow = Spow.(ctyp).pow;
    Sbsl.(ctyp).mean = mean(pow, 3);
    Sbsl.(ctyp).std = std(pow, 0, 3);
end
% Save it
save(Sdc.mat_out, 'Sbsl');
% Save info too to get param to know if computation with same param was already
% done
info = Sdc.info; %#ok
save([Sdc.dir, fsep, 'info.mat'], 'info');

% Do baseline normalization
function Spow = bsl_norm(Spow, Sbsl)
[typ, Nty] = get_names(Spow);

for i = 1 : Nty
    ctyp = typ{i};
    Sp = Spow.(ctyp);

    Nt = length(Sp.time);

    mBL = repmat(Sbsl.(ctyp).mean, [1 1 Nt]);
    sBL = repmat(Sbsl.(ctyp).std, [1 1 Nt]);

    powr = Sp.pow;
    Sp.pow = (powr - mBL) ./ sBL;

    % Add mean power across ROI for further figures
    Sp.mean_roi = pow_mean_roi(Sp);
    
    Spow.(ctyp) = Sp;
end

% Normalize and save data 
function save_pow(Spow, Sdc)

save(Sdc.mat_out, 'Spow');
info = Sdc.info; %#ok
pdir = Sdc.dir;
save ([pdir, fsep, 'info.mat'], 'info')

%- Figure of power across trials for each ROI
% only of normalized power
if ~Sdc.norm_bsl.do
    return;
end

pfig = make_dir([pdir, fsep, 'powmat_fig']);
[typ, Nty] = get_names(Spow);
for i = 1 : Nty
    ctyp = typ{i};
    cp_inv_powmat_fig(Spow.(ctyp).mean_roi, pfig)
end


% Mean pow accross ROI for figures
function Smean = pow_mean_roi(Spowz)
ins = Spowz.inside;
atlas = Spowz.atlas;
pow = Spowz.pow;
[ntr, ~, Nt] = size(pow);

% All labels for dipoles
alab = atlas.label;
% Unique ROI name
ulab = unique(alab);
ulab = ulab(~cellfun(@isempty, ulab)); 
Na = length(ulab);

Smean = Spowz;
% Mean across dipole of the same ROI
meanroi = zeros(Na, ntr, Nt);

% Keep dipole number / ROI
Ndip = zeros(Na, 1);

for i = 1 : Na
    lab = ulab{i};
    ilab = strcmp(alab, lab) & ins;
    Ndip(i) = sum(ilab);
    powr = pow(:, ilab, :);
    meanroi(i, :, :) = squeeze(mean(powr, 2));
end
Smean.pow = meanroi;
Smean.mpow = squeeze(mean(Smean.pow, 2));
Smean.label = ulab;

labf = cell(Na, 1);
for i = 1 : Na
    lab = ulab{i};
    ii = find(strcmp(alab, lab)==1, 1, 'first');
    labf{i} = [atlas.fullname{ii}, ' ', atlas.hemis{ii}];
end
Smean.labfull = labf;
Smean.Ndip = Ndip; 

% Subdirectory name according to time parameters
function tnam = dir_toi(toi, tpre)
if nargin < 2 || isempty(tpre)
    tpre = 'toi';
end
ti = num2str(toi(1), '%1.3f');
tf = num2str(toi(end), '%1.3f');
dt = num2str(toi(2)-toi(1), '%1.3f');
tnam = name_save(strrep([tpre, '_', ti, '_', tf, '_dt', dt], '.', ''));

% Prepare dics parameters and data paths
function Sdp = prep_dics(bopt)

%----- Check for cond_param parameter
cond_param = bopt.param.cond_param;
% Be sure cond_in are member of acond (avalaible conditions)
isn = ~ismember(cond_param(:, 1), bopt.acond);
if any(isn)
    if all(isn)
        warning('No condition set in cond_param is available in the cleaned dataset -- %s\n', bopt.info);
        return
    else
        warning('%d condition(s) not available in cleaned dataset -- %s', sum(isn), bopt.info);
        warning(' --> %s\n', strjoint(cond_param(isn, 1), ' ; '));
    end    
    cond_param = cond_param(isn, :);    
end

% Make cond_out names compatible with file and folder names to be crated and
% whith structure field names
cond_param(:, 2) = cellfun(@name_save, cond_param(:, 2), 'UniformOutput', 0);

%----- Condition for normalization
cbsl = bopt.param.norm_cond_out;
if ~isempty(cbsl)
    if iscell(cbsl)
        cbsl = cbsl{1};
    end
    cout = cond_param(:, 2);
    isb = ismember({cbsl}, cout);
    if ~any(isb)
        cbsl = [];
    else
        cbsl = cout{isb};
        ib = find(isb==1);        
        nc = 1 : length(cout);
        % Order condition to have baseline as first
        iord = [ib setxor(nc, ib)];
        cond_param = cond_param(iord, :);
    end 
end
cond_in = cond_param(:, 1);
cond_out = cond_param(:, 2);
cond_toi = cond_param(:, 3);

% Number of output data
Nc = length(cond_in);

%---- Frequency parameters for ft_freqanalysis
%%% TO DO: for now, FOI with only one value
ofreq = prep_freq_param(bopt.param.freq_param, bopt.dir);

% Number of frequency analysis
Nf = length(ofreq);

% Compare new preproc parameters with previous one use for power computation
preproc = loadvar([fileparts(bopt.data), filesep, 'preproc.mat']);

% New clean data have been computed just before --> need to do a new power
% computation
new_clean = bopt.new_clean;

Sdp = [];
param_bsl = [];
pbsl = [];
% Loop over frequency analysis
for i = 1 : Nf
    
    ofq = ofreq(i);    
    cfg = ofq;
    
    % Keep param/info
    param = ofq.param; 
    fqnam = param.name;
    
    isbsl = ~isempty(cbsl);
    % Loop over conditions
    for j = 1 : Nc
        cdin = cond_in{j};
        cdout = cond_out{j};       
        
        % TOI paramaters
        toi = cond_toi{j};
        cfg.toi = toi;
        param.toi = toi;      
        
        % Subdirectory depending on condition name and toi param
        pcond = make_dir([ofq.dir, filesep, cdout]);
        pres = make_dir([pcond, filesep, dir_toi(toi)]);
        
        spre = 'Spow';
        if ~isempty(cbsl)           
            % normz subdirectory 
            if j>1
                pnorm = make_dir([pres, fsep, 'normz']);
                % bsl indication directory
                pres = make_dir([pnorm, fsep, dir_toi(cond_toi{1}, 'bsl')]);  
                isbsl = 0;
            else
                spre = 'Sbsl';
            end
        end
        % Output MAT filepath
        pmat = [pres, fsep, spre, '_', cdout, '.mat'];
        
        % Associated info
        info = [];
        info.param = param;
        info.subj = bopt.info;
        info.cond_in = cdin;
        info.cond_out = cdout; 
        info.preproc = preproc;
               
        % Data for baseline computation
        if isbsl
            % Check if baseline already compute (with same parameters)
            param_bsl = info;          
        end
        info.param_bsl = param_bsl;
        
        % Check if computation has been already done     
        isdone = is_done(pres, info, new_clean);
        
        if isdone 
            % Keep precompute baseline
            if isbsl
                pbsl = pmat;
            end
            continue
        end
        
        % Store param & info to do dics computation
        Sp = [];
        Sp.cond_in = cdin;
        Sp.cfg = cfg;
        Sp.isbsl = isbsl;
        
        Sp.info = info;
        Sp.mat_out = pmat;
        Sp.dir = pres;
        Sp.norm_bsl.do = ~isempty(cbsl);
        % Store precomputed baseline
        Sp.norm_bsl.precomp = loadvar(pbsl);    
        
        Sdp.(fqnam).(cdout) = Sp;
    end
end

% Check if computation have been previously done with exactly the same
% parameters
function isdone = is_done(pmat, info, new_clean)

isdone = 0;

if new_clean
    return;
end

pdir = fileparts(pmat);

% Redo Spow computation if any parameters have been changed
ppr = [pdir, fsep, 'param_info.mat'];
if ~exist(ppr, 'file') ||...
        ~isequal(info, loadvar(ppr)) ||...
        ~exist(pmat, 'file')
    return;
end

% Spow or Sbsl exists with required parameters
isdone = 1;


% Run DICS only if required parameters are set
function isok = check_bopt(bopt)
isok = ~any(cellfun(@isempty, {bopt.param.freq_param,...
                        bopt.param.cond_param,...
                        bopt.acond,...
                        bopt.fwd,...
                        bopt.data}));