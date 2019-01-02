function  allTrials = extract_trials(ftData, Spar, eopt)
% Extract trial from continuous dataset ftData
%-CREx180726

% Trial number association by condition
isa_t = strcmp(eopt.rm_trials_cond, 'same');
% eopt.datafile is required in cmeg_extract_trials by ft_definetrial if the default trialfun
% function is used ('ft_trialfun_general')
eopt.datafile = Spar.dir.raw;
% Return event structures without trials falling inside strong artefacts (if any)
[Strig, allev, iwart] = prepare_event(Spar, eopt.trigfun);

cond = {Strig(:).name}';
Nc = length(cond);

for k = 1 : Nc
    cnam = cond{k};
    eopt.trig = Strig(k);   
    % Associated pre and post stimulation times for epoching
    eopt.prestim = eopt.trig.dt_s(1);
    eopt.poststim = eopt.trig.dt_s(2);
    trials = cmeg_extract_trials(ftData, allev, eopt);

    % Identify trials that contain portion of artefact
    isart = artefact_inside(trials.sampleinfo, iwart);
    if isa_t
        if k==1
            allart = zeros(length(isart), Nc);
        end
        allart(:, k) = isart;
    end
    % Be sure to set grad in mm for source analysis
    trials.grad = ft_convert_units(trials.grad, 'mm');
    %-- Resampling   
    if ~isempty(eopt.res_fs)
        cfg = [];
        cfg.resamplefs = eopt.res_fs;
        trials = ft_resampledata(cfg, trials); 
    end
    trials.hdr.artefact = isart;
    allTrials.(cnam) = trials;
end
% Same artefact info for all conditions
if isa_t && any(allart(:))
    allart = sum(allart, 2) > 0;
    for k = 1 : Nc
        allTrials.(cond{k}).hdr.artefact = allart;
    end
end

function [Strig, allev, iart] = prepare_event(Spar, trigfun)

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
    iart = [];
    return;
end
% Artefact windows
Fs = hdr.Fs;
iart = wina.*Fs;

% Identify trials which contain a portion of padded artefact
function isart = artefact_inside(iwtr, iwart)
Ntr = length(iwtr(:, 1));
isart = zeros(Ntr, 1);
if isempty(iwart)
    return;
end

Na = length(iwart(:, 1));
for i = 1 : Ntr
    intr = iwtr(i, 1) : iwtr(i, 2);
    for j = 1 : Na
        inta = iwart(j, 1) : iwart(j, 2);
        if any(ismember(intr, inta))
            isart(i) = 1;
            continue;
        end
    end
end

fprintf('\nKept trials outside strong artefacts: %d/%d\n', sum(~isart), Ntr);