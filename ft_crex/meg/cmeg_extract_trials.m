function trials = cmeg_extract_trials(ftData, Sev, opt)
% Extract trials with ft_definetrial and ft_redefinetrial 
%
%-- ftData** : the data structure asoutput by feildtrip functions
%
%-- Sev** : the events structure (output of ft_read_event)
%
%-- opt : option structure with fields:
%
%  - opt.trialfun : name of specific function to keep or reject trials 
%       [default :'ft_trialfun_general']
%
%  - opt.datafile(**) : data path, mandatory if the default trialfun is
%  used
%    
%  - opt.prestim : time before stimulus (s) [default : 0.5]
% 
%  - opt.poststim : time after stimulus (s) [default : 1]
%
%  - opt.trig : structure of trigger definition with fields :
%           value**     : trigger values 
%           eventyp   : type for events in event stucture [default : 'TRIGGER']
%           resptyp   : type for responses [default : 'RESPONSE']
%           rightresp : code of trigger for right responses if needed
%   This parameters are eventually used by the project's specific function
%   defined by opt.trialfun.
%
% (** Mandatory variables)
% ---

%---
% Check dft structure - put default value if missing field 

dopt = struct('prestim', 0.5,...
                'poststim', 1,...
                'trialfun','ft_trialfun_general',...
                'trig', struct('eventyp','TRIGGER','resptyp','RESPONSE'),...
                'datafile', []);

opt = check_opt(opt, dopt); 

% fsample needed for custom determination of trials to keep
if isfield(ftData,'fsample')
    fsample = ftData.fsample;
elseif (isfield(ftData,'hdr') && isfield(ftData.hdr,'Fs'))
    fsample = ftData.hdr.Fs;
else
    disp('fsample not found...')
    tim = ftData.time{1};
    fsample = (length(tim)-1)./(tim(end) - tim(1));
end

cfg = [];
cfg.event = Sev;
cfg.trialfun = opt.trialfun; 
% If trialfun is the Fieldtrip default one 'ft_trialfun_general', the
% datafile path is required to extract events (l. 74)
if (strcmp(opt.trialfun, 'ft_trialfun_general') == 1)
    if ~isempty(opt.datafile)
        cfg.headerfile = opt.datafile;
    else
        warning(['Datafile path required to run extraction with default',...
        'ft_trialfun_general function'])
    end
end
% Trigger codes
cfg.trialdef.eventvalue = opt.trig.value; 
cfg.trialdef.prestim    = opt.prestim;
cfg.trialdef.poststim   = opt.poststim;
% Name of the trigger channel to extract
cfg.trialdef.eventtype  = opt.trig.eventyp;  
cfg.trialdef.resptype  = opt.trig.resptyp;

if isfield(opt.trig,'rightresp') && ~isempty(opt.trig.rightresp)
    cfg.trialdef.rightresp  = opt.trig.rightresp;
end

cfg.fsample = fsample;

cfg_trial   = ft_definetrial(cfg);

% Ajout verification si longueur donnees > indice des essais dans
trl = cfg_trial.trl;

if isfield(opt.trig, 'tshift') && ~isempty(opt.trig.tshift) &&...
        opt.trig.tshift~=0
    tshi = opt.trig.tshift;
    fprintf('\nTrial shift : %1.3f\n', tshi);
    Nshi = round(tshi*fsample);
    trl(:, 1:2) = trl(:, 1:2) + Nshi;
end
    
Ns = length(ftData.time{1});
if any(trl(:, 2) > Ns)
    fprintf('\n\n!!! Continuous data length inferior to trials indices definition')
    fprintf('Data length : %d samples\n', Ns)
    fprintf('Trials that fall outside:\nOnset    End\n')
    disp(num2str(trl(trl(:,2)>Ns, 1:2)))
    disp('---')
    Nini = length(trl(:,1));
    trl_cr = trl(trl(:,2) < Ns, :);
    fprintf('Keeping %d / %d trials\n', length(trl_cr(:,1)), Nini);
    trl = trl_cr;
end

cfg = [];
cfg.trl = trl;  
trials = ft_redefinetrial(cfg, ftData); 

