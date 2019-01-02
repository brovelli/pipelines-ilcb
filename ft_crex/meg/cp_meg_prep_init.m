function Sdb = cp_meg_prep_init(Sdb, opt)
% Initialize preprocessing parameters
%- path to previous txt files holding preprocessing parameters 
% in db_fieldtrip/PROJ/(group)/SUBJ/meg/preproc_param directory (bad sensors, strong artefacts,
% bad trials)
%- check if cleanup fig are to be processed (new MEG data)
%- check if ICA is to be processed (new sensor or strong artefact removing)
%
%-CREx-180726

% Init param_txt paths -> Sdb.meg.preproc.param_txt
Sdb = init_param_txt(Sdb, opt);

% Add the previously selected BAD elements read from TXT files
% (sensors, strong artefact windows, ICA components, trials)
% -> Sdb.meg.preproc.param_run
% + check if cleanup_fig and/or ICA is to be done -> Sdb.meg.preproc.do
Sdb = init_preproc(Sdb, opt);

% Initialize all default paths for param_txt files
function Sdb = init_param_txt(Sdb, opt)

Ns = length(Sdb);
for i = 1 : Ns
    Smeg = Sdb(i).meg;
    dp_par = Smeg.txt_dir;
    dpp = [dp_par, fsep]; 
    
    dpinfo = Smeg.info;
  
    drun = Smeg.run.dir;
    Nr = length(drun);
    
    Sprun = cell(Nr, 1);
    Spar = []; 
    Spar.dir = dp_par;
    % Check for rmsens_allrun.txt file
    Spar.rms.allrun = [dpp, 'rmsens_allrun.txt'];
    for j = 1 : Nr
        srun = drun{j};
        % RM sensors / run
        Spar.rms.(srun) = [dpp, 'rmsens_', srun, '.txt'];
        % RM strong artefact / run
        Spar.rma.(srun) = [dpp, 'rmart_', srun, '.txt'];
        
        % RM ICA component / run
        Spar.rmc.(srun) = [dpp, 'rmcomp_', srun, '.txt'];
        
        % RM trials depending on conditions
        
        % Prepare markers name + values + (type) according to hdr_event.event
        % values and trigfun function
        cond = get_cond(dpinfo{j}, opt);

        Nc = length(cond);
        
        % Same trials for conditions == small shifting // same trigger
        Spar.rmt.(srun).allcond = [dpp, 'rmtrials_', srun, '_allcond.txt'];
        % RM trials / each conditions
        for k = 1 : Nc
            scond = cond{k};
            Spar.rmt.(srun).(scond) = [dpp, 'rmtrials_', srun, '_', scond, '.txt'];
        end
        % Condition to extract
        Sprun{j}.conditions = cond;
        % Associated epoching time [prestim postim]
        Sprun{j}.dt_s = epoching_opt(opt, cond);
    end
    Sdb(i).meg.preproc.param_txt = Spar;
	Sdb(i).meg.preproc.param_run = Sprun;
end
 
%- Define all conditions depending on what is returned by the trigfun function
% This function may need sometimes the event structure to set conditions
% depending on the markers recording during the run
function cond = get_cond(pinfo, opt)
cond = opt.epoched.conditions;
if ~isempty(cond)
    return;
end
hdr = loadvar([pinfo, fsep, 'hdr_event']);
ev = hdr.event;
ftrig = str2func(opt.epoched.trigfun); 
Strig = ftrig(ev);
cond = {Strig(:).name};

% Set dt_s consistant with conditions number (ex. if dt_s is the same for
% all contidions, only one vector [1x2] has been set as input parameter)
function dt_s = epoching_opt(eopt, cond)
Nc = length(cond);
dt_s = eopt.epoched.dt_s;
if Nc > 1 && length(dt_s(:, 1))==1
    dt_s = repmat(dt_s, Nc, 1);
end
% Be sure all times > 0 s
dt_s = abs(dt_s);

% Check if data visualization and ICA was already done and read the param in TXT
% files
function Sdb = init_preproc(Sdb, opt)

isa_s = strcmp(opt.continuous.rm_sens_run, 'same');
isa_t =  strcmp(opt.epoched.rm_trials_cond, 'same');

fica =  opt.continuous.ica.reject;

% Init data process directory
sfilt = meg_prep_dir(opt);

% Loop on subject to set preprocessing parameters
Ns = length(Sdb);
for i = 1 : Ns
    
    % Subject ID
    idsubj = Sdb(i).sinfo;
    
    % MEG structure
    Smeg = Sdb(i).meg;
    
    % Raw data paths
    dp_raw = Smeg.raw;
    
    % Preprocessing directory for continuous data
    dp_prep = Smeg.preproc.dir;
    
    % Clean data directory <--> where the clean data set will be saved for analysis
    dp_clean = Smeg.analysis.dir;
    
    % info folder when hdr_event have been stored
    dp_info = Smeg.info;
    
    % Folder with parameter TXT files
    partxt = Smeg.preproc.param_txt;
    
    % Run directory
    rdir = Smeg.run.dir;
    Nr = length(rdir);
    
    % Preprocessing parameters for each run
    Srun = Smeg.preproc.param_run;
    
    % Sdo will hold the flag to indicate which new preprocessing computation
    % have to be done
    Sdo = [];
    % New bad sensor selection or new cleanup fig 
    Sdo.rms = zeros(Nr, 1);
    % New artefact rejection => new dispdata and fftplots are
    % to do in cleanup_fig directory
    Sdo.rma = zeros(Nr, 1);
    % New ICA if fica==1
    Sdo.ica = zeros(Nr, 1);
    % Bad component selection
    Sdo.rmc = zeros(Nr, 1);
    % New bad trial interactive identification (by ft_rejectvisual)
    Sdo.rmt = zeros(Nr, 1);
    % Prepare the epoched data for source analysis with the real opt
    % parameters for filtering and resampling + the RM sensors / comp /
    % trials 
    Sdo.epch = zeros(Nr, 1);
    
    % Loop over runs
    for j = 1 : Nr
        % Update run statut (that could have changed during the final review)
        Smeg.run.valid(j) = update_valrun(Smeg.run.valtxt{j});

        % If run is not valid, keep all do flag = zeros
        if ~Smeg.run.valid(j)
            continue;
        end
        
        % Run directory ("run_*")
        srun = rdir{j};
        sinfo = [idsubj, ' -- ', srun];
        
        % Set paths to meg directories
        Spar = Srun{j};
        pprep = dp_prep{j};
        Spar.dir = struct('raw', dp_raw{j},...
                            'info', dp_info{j},...
                            'preproc', pprep,...
                            'cleanup_fig', make_dir([pprep, fsep, 'cleanup_fig']),...
                            'ica', [pprep, fsep, 'ica'],...
                            'clean', make_dir([dp_clean{j}, fsep, sfilt]));
         
        % Default folder for preprocessing figures
        pfig_data = [Spar.dir.cleanup_fig, fsep, 'datadisp'];
        pfig_fft = [Spar.dir.cleanup_fig, fsep, 'fftplots'];       
           
        % Figure with superimposed fft for bad channel manual selection
        Spar.rms_fig = [pfig_fft, fsep, 'fftstack_interact.fig'];    
        
        %--- Read all previously saved parameters in txt files according to
        % paramaters isa_s, fica and isa_t
        Spar = read_param_txt(Spar, partxt, srun, isa_s, fica, isa_t);

        %------------------ Do STRONG ARTEFACT detection
        % Do artefact detection for the first time (except if it was done previously)
        if ~exist(partxt.rma.(srun), 'file')            
            Sdo.rma(j) = 1;            
        end           
        %------------------ Do BAD SENSOR identification
        % Do data figures and bad sensor identification for the first time
        if ~exist(pfig_data, 'dir') || ~exist(pfig_fft, 'dir')...
                || ~exist(Spar.rms_fig, 'file') || ~exist(partxt.rms.(srun), 'file')          
            Sdo.rms(j) = 1;            
        end        
        
        %------------------ Do ICA 
        if fica 
            pica = make_dir(Spar.dir.ica);
            if (~exist([pica, fsep, 'icaComp_res.mat'], 'file') ||...
                ~exist([pica, fsep, 'preproc_ica.mat'], 'file'))
                % Need to make ICA and/or redefine preproc_ica.mat
                Sdo.ica(j) = 1;
                Sdo.rmc(j) = 1;
            else
               [Sdo.ica(j), Spar] = check_ica_preproc(Spar, partxt, srun, sinfo); 
               if Sdo.ica(j)
                   Sdo.rmc(j) = 1;
               end
            end
            % Ask for bad component with already computed ICA
            if ~Sdo.ica(j) && ~exist(partxt.rmc.(srun), 'file')
                Sdo.rmc(j) = 1;
                % Figures are required to select bad components
                pica_fig = [pica, fsep, 'ica_fig'];
                if ~exist(pica_fig, 'dir') ||...
                        isempty(dir([pica_fig, fsep, '*.png']))
                    Sdo.ica(j) = 1;
                end
            end            
        end
        
        %------------------- Do BAD TRIALS identification
        % Check if bad trials identification have already be
        % done previously (depending on rms, rmc, rma and opt.epoched
        % parameters (should be the same for dt and conditions name)
        if Sdo.rms(j) || Sdo.ica(j) 
            Sdo.rmt(j) = 1;
        else
            [Sdo.rmt(j), Spar] = check_prep_trials(Spar, isa_t);
        end
        
        %------------------- Do the final epoching !!
        if any([Sdo.rms(j) Sdo.rma(j) Sdo.ica(j) Sdo.rmt(j)])
            Sdo.epch(j) = 1;
        else
            Sdo.epch(j) = check_epochs(Spar, isa_t);
        end
        
        % Add the cleanTrials.mat path to the Sdb
        if ~Sdo.epch(j)
            Smeg.analysis.clean_mat{j} = [Spar.dir.clean, fsep, 'cleanTrials.mat'];           
        else
            Smeg.analysis.new_clean(j) = 1;
        end

        Smeg.analysis.clean_dir{j} = Spar.dir.clean;
        
        Smeg.preproc.param_run{j} = Spar;    
    end
    Smeg.preproc.do = Sdo;
    Sdb(i).meg = Smeg;
end

% Check for preprocessing parameters applied to previously saved epoched data
function isep = check_epochs(Spar, isa_t)
isep = 0;
pep = Spar.dir.clean;
ppr = [pep, fsep, 'preproc.mat'];
pdat = [pep, fsep, 'cleanTrials.mat'];

if ~exist(ppr, 'file') || ~exist(pdat, 'file')
    isep = 1;
    return;
end

prep = loadvar(ppr);
prep = prep.rm;
Srm = Spar.rm;
isame = [are_same(Srm.sens, prep.sens)
    are_same(Srm.comp, prep.comp)
    are_same(Srm.art, prep.art)];
if ~all(isame)
    isep = 1;
    return;
end

if isa_t
    isame = are_same(Srm.trials.allcond, prep.trials.allcond);
    if ~isame
        isep = 1;
    end
    return;
end

cond = Spar.conditions;
Nc = length(cond);
for k = 1 : Nc
    cnam = cond{k};
    if ~isfield(prep.trials, cnam)
        isep = 1;
        break;
    end
    isame = are_same(Srm.trials.(cnam), prep.trials.(cnam));
    if ~isame
        isep = 1;
        break;
    end
end

function [isrmt, Spar] = check_prep_trials(Spar, isa_t)
isrmt = 0;
pptr = [Spar.dir.preproc, fsep, 'preproc_trials.mat'];
if ~exist(pptr, 'file')
    isrmt = 1;
    return;
end
% Check if preprocing parameters were the same when trials were identified
prept = loadvar(pptr);

Srm = Spar.rm;
isame = [are_same(Srm.sens, prept.rm.sens)
    are_same(Srm.comp, prept.rm.comp)
    are_same(Srm.art, prept.rm.art)];
if ~all(isame)
    isrmt = 1;
    return;
end
% Check for identical conditions (at least the one required
% are inside prept.condition) and dt (should be the same !
cond = Spar.conditions;
Nc = length(cond);
for k = 1 : Nc
    cnam = cond{k};
    ip = strcmp(cnam, prept.conditions);
    if ~sum(ip)
        isrmt = 1;
        break;
    else
        dt_s = Spar.dt_s(k, :);
        dt_sp = prept.dt_s(ip, :);
        if any(dt_s - dt_sp) > 0
            isrmt = 1;
            break;
        end
    end
end
if isa_t 
    Spar.Ntr.allcond = prept.Ntr.(cnam);
else
    Spar.Ntr = prept.Ntr;
end
function iss = are_same(par1, par2)
isi = sum([isempty(par1) isempty(par2)]);
if isi==2
    iss = 1;
    return;
elseif isi==1
    iss = 0;
    return;
end
par1 = prep_param(par1);
par2 = prep_param(par2);
iss = isempty(setxor(par1, par2));

% Read all param previously written in txt files
function Spar = read_param_txt(Spar, partxt, srun, isa_s, fica, isa_t)
      
%--- RM Sensors
if isa_s 
    frun = 'allrun';
else
    frun = srun;
end
% Same sensors selection to delete for all runs
Spar.rm.sens = read_param(partxt.rms.(frun), 'sens');  

%--- RM Artefacts and reject Component if ICA
if fica
    % Read artefacts
    Spar.rm.art = read_param(partxt.rma.(srun), 'art');
    % Read components
    Spar.rm.comp = read_param(partxt.rmc.(srun), 'comp');
else
    Spar.rm.art = [];
    Spar.rm.comp = [];
end

%--- RM trials
% Same trials for all conditions (== conditions designed the same trial
% with slighlty different epoching interval - cf. Stim/Action/Reward
% design type)                  
cond = Spar.conditions;
if isa_t
    fn = {'allcond'};
else
    fn = cond;
end
Nc = length(fn);
for k = 1 : Nc
    cnam = fn{k};
    ptr = partxt.rmt.(srun).(cnam);
    Spar.rm.trials.(cnam) = read_param(ptr, 'trl');
end
        
% Read BAD element selection in TXT file
function parm = read_param(pfile, type)
if ~exist(pfile, 'file')
    parm = [];
    return;
end
fid = fopen(pfile);
switch type
    case 'sens'
        sp = textscan(fid, '%s');
        parm = sp{1};
    case 'art'
        sp = textscan(fid, '%f\t%f');
        parm = [sp{1} sp{2}];
    otherwise
        sp = textscan(fid, '%d');
        parm = sp{1};
end
fclose(fid);      

% Check for previous ICA use regarding to new preprocessing options (rm_sens and
% rm_art)
function [isnew, Spar] = check_ica_preproc(Spar, partxt, srun, sinfo)
   
%-----------
% Check if we can use the previous ICA
% ==> requires the same bad channels and strong artefact
% selection
%-----------
rms = Spar.rm.sens;
rma = Spar.rm.art;
pica = Spar.dir.ica;
comp = loadvar([pica, fsep, 'icaComp_res.mat']);  
% Number of components for already processed ICA
Spar.Ncomp = length(comp.topolabel);    

ppica = [pica, fsep, 'preproc_ica.mat'];

rm_prev = loadvar(ppica); 
% Check for parameters change // those required to keep ICA
rms_mod = check_mod(rms, rm_prev.rms);
% Check if new rma has been added / removed in param_txt file
% compare to initial rma applied when ICA was computed
rma_mod = check_mod(rma, rm_prev.rma, 'art');                

%------ Warning message for removed sensors / artefacts after ICA
isnew = warning_ica(rms_mod.removed, rma_mod.removed, sinfo, 'rem');

if isnew
    return;
end

% Add removed sensors/artefacts to the analysis
Spar.rm.sens = [rms ; rms_mod.removed];
Spar.rm.art = [rma ; rma_mod.removed];

% Write txt files
write_bad(partxt.rms.(srun), Spar.rm.sens, 'sens');
write_bad(partxt.rma.(srun), Spar.rm.art, 'art');

%------ Warning for added sensors / artefacts after ICA
isnew = warning_ica(rms_mod.new_add, rma_mod.new_add, sinfo, 'add');
if isnew
    return;
end
% Add the new sensors/artefact in preproc_ica.rms.after
% to not display ICA warning for further processings
preproc_ica = rm_prev;

preproc_ica.rms.after = [rm_prev.rms.after; rms_mod.new_add];
preproc_ica.rma.after = [rm_prev.rma.after; rma_mod.new_add];
save(ppica, 'preproc_ica');

% Check for modification of the BAD channels or artefact window selection compare to the one
% use before ICA
%-- rms_mod.new_add : list of new BAD sensors/windows manually added in the rmsens.TXT file
%   compare to the previous ones
%-- rms_mod.removed : the BAD sensors/windows that have been deleted from the TXT file
%   while there where initially used to prepare data for ICA
function rm_mod = check_mod(cpar, prm_ica, cart)
if nargin < 3
    isart = 0;
else
    isart = strcmp(cart, 'art');
end
ppar = [prm_ica.before ; prm_ica.after];
ppar = prep_param(ppar);
cpar = prep_param(cpar);

% If any difference between elements in TXT file and what was RM before or after
% ICA
rm_mod = [];
rm_mod.new_add = [];
% Removed artefact or sensor compare to the BEFORE ICA selection
% (need to recompute ICA if a channel is finally put back to the data set)
rm_mod.removed = [];

% Not the same preproc param
if ~isempty(setxor(cpar, ppar))
    rm_mod.new_add = out_param(setdiff(cpar, ppar));    
    rm_mod.removed = out_param(setdiff(prep_param(prm_ica.before), cpar));
    if isart
        if isempty(rm_mod.new_add)
            rm_mod.new_add = [];
        end
        if isempty(rm_mod.removed)
            rm_mod.removed = [];
        end
    end
end

% Format the output param (for artefact windows)
function cpar = out_param(cpar)
if ~isempty(cpar) && ~isempty(strfind(cpar{1}, '-')) 
    Na = length(cpar);
    spar = cellfun(@(x) textscan(x, '%f - %f'), cpar, 'UniformOutput', 0);
    cpar = zeros(Na, 2);
    for i = 1 : Na
        cpar(i, :) = [spar{i}{1} spar{i}{2}];
    end
end

% Prepare param to be compare (num -> cell of strings) (cf. for artefact
% windows)
function cpar = prep_param(cpar)
if ~isempty(cpar) && ~iscell(cpar)
    sz = size(cpar);
    if sz(1)==1 && sz(2) > 2
        cpar = cpar';
    end
    cpar = cellstr(num2str(cpar));
    cpar = cellfun(@(x) strjoint(strsplitt(x, ' '), ' - '), cpar, 'UniformOutput', 0);        
end

% Warning message for ICA if preprocessing parameters change is detected
% with the choice to launch a new ICA or to remove new parameters after ICA
% rejection / or to reintroduce parameters as defined when doing the ICA
function new_ica = warning_ica(crms, crma, sinfo, swar)
if ~any([~isempty(crms) ~isempty(crma)])
    new_ica = 0;
    return;
end
msg_s = warn_msg(crms, 'sensor', swar);
msg_a = warn_msg(crma, 'artefact', swar);
    
% Initialize a msgbox
smsg = [{sinfo}; msg_s ; msg_a];
if strcmp(swar, 'rem')
    new_ica = ica_msg_box('Title', 'Parameter change detected', 'string', smsg);
else
    new_ica = ica_msg_box_add('Title', 'Parameter change detected', 'string', smsg);
end

function msg = warn_msg(clist, typ, swar)
if isempty(clist)
    msg = [];
    return;
end
if ~iscell(clist)
    clist = prep_param(clist);
    clist = cellfun(@(x) ['[',x, ']'], clist, 'UniformOutput', 0);
end
if strcmp(swar, 'rem')
    msg = {[typ, 's previously removed prior to ICA calculation are no longer included'];
        ['in the list of ', typ, 's to be removed : ']; strjoint(clist, '; ')};
else
    msg = {['New ', typ, '(s) to remove from data have been added'];
    ['in the list of ', typ, 's to be removed : ']; strjoint(clist, '; ')};
end
msg{1}(1) = upper(msg{1}(1));
