function [trl, event] = trialfun_te(cfg)
% Specific function to define the trials to be kept for MEG_TE project
% depending on the triggers values returned by trigfun_te
% Only the trials that are part of a complete sequence 's' 'a' 'r' with a correct answer 
% are kept 
% Correct answer is defined as the 2nd trigger value in trigfun_te (Strig.value(2)) for
% Strig.name == 'R' 
%
%-CREx180724

% Sample frequency
Fs   = cfg.fsample;

% All events
event = cfg.event;

% Find OK events
[event, tval] = trial_sar(event, cfg.trialdef.eventvalue);
samples = [event.sample]';

% Determine the number of samples before and after the event
Npre  = -round(cfg.trialdef.prestim  * Fs);
Npost =  round(cfg.trialdef.poststim * Fs);

Ns = length(samples);
% trl for trial definition in ft_redefinetrial
trl = [samples+Npre samples+Npost repmat(Npre, Ns, 1) tval];

% Valid triggers = the one which are in a 's' 'a' 'r' complete sequence + with a
% correct response (nor late nor incorrect)
function [event, tval] = trial_sar(event, trval)

Strig = trigfun_te;

% Find the associated condition name
cval = {Strig(:).value};
icond = cellfun(@(x) all(ismember(trval, x)), cval);
tfun = Strig(icond).fun;

% Only keep trial with a complete sequence 'S' 'A' 'R'
cnam = {Strig(:).name}';
Strig = Strig([find(strcmp(cnam, 'S')) find(strcmp(cnam, 'A')) find(strcmp(cnam, 'R'))]);

val   = [event.value]';
Nev  = length(event);
      
% Stimulus (nstim)
Ne = length(Strig);

cind = cell(3, 1);
for k = 1 : Ne
    mval = Strig(k).value;
    n_mk = length(mval);
    mind = find( sum((repmat(val, 1, n_mk)==repmat(mval, Nev, 1)), 2) == 1 );
    cind{k} = [mind ones(length(mind), 1)*k];
end

% Group indices and sort
ind  = cell2mat(cind);
[ind(:, 1), ia]  = sort(ind(:,1));
ind(:, 2) = ind(ia, 2);

%-- Select sequences of stim, act, rew = 1,2,3 in ind(:,2)
iok = strfind(ind(:, 2)', [1 2 3])';

ind_ok = [ind(iok, 1) ind(iok+1, 1) ind(iok+2, 1)]; 

%-- Only keep Correct trial, when Response == 722 (== Strig(3).value(2) as
%   defined in trigfun_te for trigger of type 'R'
irok = val(ind_ok(:, 3)) == Strig(3).value(2);

fprintf('Kept trials: %d/%d\n', sum(irok), length(irok));
ind_ok = ind_ok(irok, :);

ind_ok = ind_ok(:);
%-- Now, keep trigger with trval value
Ni = length(ind_ok);
Nv = length(trval);
ind_tr = ind_ok(sum((repmat(val(ind_ok), 1, Nv)==repmat(trval, Ni, 1)), 2) == 1);

event = event(ind_tr);

% At least, keep trigger value column (but with value that is transfomed by
% Strig(icond).fun anonymous function
tval = tfun([event.value])';
