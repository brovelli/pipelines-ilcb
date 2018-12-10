function  allTrials = extract_trials(ftData, Spar, eopt)
% Extract trial from continuous dataset ftData
%-CREx180726

% eopt.datafile is required in cmeg_extract_trials by ft_definetrial if the default trialfun
% function is used ('ft_trialfun_general')
eopt.datafile = Spar.dir.raw;
% Return event structures without trials falling inside strong artefacts (if any)
[Strig, allev] = event_padart(Spar, eopt.trigfun);

cond = {Strig(:).name}';
Nc = length(cond);
for k = 1 : Nc
    cnam = cond{k};
    eopt.trig = Strig(k);   
    % Associated pre and post stimulation times for epoching
    eopt.prestim = eopt.trig.dt_s(1);
    eopt.poststim = eopt.trig.dt_s(2);
    trials = cmeg_extract_trials(ftData, allev, eopt);

    % Be sure to set grad in mm for source analysis
    trials.grad = ft_convert_units(trials.grad, 'mm');
    %-- Resampling   
    if ~isempty(eopt.res_fs)
        cfg = [];
        cfg.resamplefs = eopt.res_fs;
        trials = ft_resampledata(cfg, trials); 
    end
    allTrials.(cnam) = trials;
end

function [Strig, allev] = event_padart(Spar, trigfun)

cond = Spar.conditions;

% Get the events structure
hdr = loadvar([Spar.dir.info, filesep, 'hdr_event']);
allev = hdr.event;
ftrig = str2func(trigfun); 
Strig = ftrig(allev);
allcond = {Strig(:).name}';
% Keep the same order as in cond
[~, ia] = ismember(cond, allcond);
Strig = Strig(ia(ia~=0));
% Add epoching time info / condition
Nc = length(Strig);
for k = 1 : Nc
    Strig(k).dt_s = Spar.dt_s(k, :);
end
% Remove artefact
wina = Spar.rm.art;
if isempty(wina)
    return;
end
% Add additional 1 s around the padart window to remove events in addition
% to them that fall inside the artefact window itself
Fs = hdr.Fs;
wina = wina.*Fs;
winpad = [wina(:,1)-Fs wina(:,2)+Fs];

% Remove all events that fall inside the artefact enlarged windows
allev = set_padart_event(allev, winpad);


% Change event type (TRIGGER or RESPONSE) to INSIDE_PADART
% in order to exclude these events for further processing (epoching)
function ev = set_padart_event(ev, winpad)

sval = cell2mat({ev.sample})';

Na = length(winpad(:,1));
for k = 1 : Na
    ide = find(sval >= winpad(k,1) & sval <= winpad(k,2));
    if ~isempty(ide)
        for j = 1 : length(ide)
            ev(ide(j)).type = 'INSIDE_PADART';
        end
    end
end
ikeep = strcmp({ev.type}', 'INSIDE_PADART') == 0;
Nini = length(ev);
ev = ev(ikeep);

fprintf('\nKept trials outside strong artefacts: %d/%d\n', sum(ikeep), Nini);