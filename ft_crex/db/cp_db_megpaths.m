function Smeg = cp_db_megpaths(praw_meg, pmeg, crun)
%-- Prepare MEG data paths from the meg directory in FIELDTRIP_DATABASE
% The MEG data directories as prepared by cp_db_import should have this
% structure:
% meg /  
%   --- (run_*)
%               --- _info
%               --- _preproc
%               ---  analysis
%
% The _info folder must contain the raw_relpath.mat file that holds the relative
% path to the subject's raw data directory in the MEG_DATABASE. This file is
% defined during cp_db_import execution.
%
% The hdr_event structure is created and saved in each _info directory as a hdr_event.mat file.
% This structure holds general information from the raw data as returned by ft_read_header, 
% the events read by ft_read_event and the headshape if a headshape file was found.
% 
% praw_meg: path to raw MEG data in DB_MEG
% pmeg: path to MEG results in DB_FT
% drun: list of run to process (as run directory name(s) initially given in
% Sdir.meg_run parameter to find raw data run directories)
%-CREx-181106

Smeg = [];
if ~exist(praw_meg, 'dir')
    return
end

if nargin < 3
    crun = [];
else
    if ischar(crun)
        crun = {crun};
    end
end

% Run directories previously imported in pmeg (excluding run_concat)
drun = rundir(pmeg, crun);

%%% TO DO : consider the run_concat if concatenation across run is required
if isempty(drun)
    warning('No MEG run folder found in %s\n', pmeg)
    return
end
% Preprocessing parameters folder
ppar = make_dir([pmeg, fsep, '_preproc_param']);

% Number of run
Nd = length(drun);
%--- Initialize data paths
%_ Directory with the preprocessing parameter txt files
Smeg.txt_dir = ppar;
%_Run directories list
Smeg.run.dir = drun;
%_ Indication if run dataset is valid for analysis
Smeg.run.valid = zeros(Nd, 1);
%_ Associated txt file containing validity flag (could be change by hand in
%_preproc_param directory)
Smeg.run.valtxt = cell(Nd, 1);
%_Number of runs
Smeg.Nrun = Nd;
%_Raw data path in DB_MEG
Smeg.raw = cell(Nd, 1);
%_Info directory holding relatave raw data path and hdr_event.mat for each run
Smeg.info = cell(Nd, 1);
%_ hdr_event data
Smeg.hdr_event = cell(Nd, 1);

%_Data preprocessing info
%  -directory ('_preproc') for each run
Smeg.preproc.dir = cell(Nd, 1);
%_ MEG analysis info
%  -main directory ('analysis')
Smeg.analysis.dir = cell(Nd, 1);
%  -clean dataset directory 
Smeg.analysis.clean_dir = cell(Nd, 1);
%  - clean data file paths for source analysis
Smeg.analysis.clean_mat = cell(Nd, 1);
%  - indication if clean data are new (= source analysis has to be done even if
%   source results already exists (Spow mat file)
Smeg.analysis.new_clean = zeros(Nd, 1);
for j = 1 : Nd
    rnam = drun{j};
    % Run directory in subject's meg folder in Fiedltrip database (previously
    % created by cp_db_import)
    prun = [pmeg, filesep, rnam];
    % _info directory that was created by cp_db_import
    pinfo = [prun, filesep, '_info'];
    % Load relative raw dataset path from db_meg subject's folder
    prel = loadvar([pinfo, filesep, 'raw_relpath.mat']);
    % Full raw data path
    if ~isempty(prel)
        praw = [praw_meg, filesep, prel];
    else
        praw = praw_meg;
    end

    % Run information (validity)
    pval = [ppar, fsep, 'validity_', rnam, '.txt'];
    Smeg.run.valtxt{j} = pval;
    Smeg.run.valid(j) = add_isval(pval);
    
    % Add hdr info in _info directory
    Smeg.info{j} = pinfo;
    Smeg.hdr_event{j} = add_hdr(praw, pinfo);    
    
    % Add run directories data paths
    Smeg.raw{j} = praw;
    Smeg.preproc.dir{j} = make_dir([prun, filesep, '_preproc']);
    % Add analysis directory
    Smeg.analysis.dir{j} = make_dir(fullfile(prun, 'analysis'));    
end

% To see later if necessary (depending on the data processing pipeline)
% Smeg.epoched.concat = make_dir(fullfile(pmeg, 'epoched', 'run_concat'));

% Find run directories (except for the 'run_concat' directory)
function drun = rundir(pmeg, crun)
if isempty(crun) ||...
        (iscell(crun) && any(cellfun(@(x) strcmp(x, '*'), crun)))
    crun = {'*'};
end
% List of run directories (excepted 'run_concat')
dlist = dirfold([pmeg, filesep, 'run_*'], {'run_concat'});
if strcmp(crun{1}, '*')
    drun = dlist;
    return
end
% Only choose specific run directories (matching with crun= cell of matching
% string with the run directory names to find in the DB_MEG database
Ns = length(crun);
drun_all = cell(Ns, 1);
% Available runs list
davr = cellfun(@(x) x(5:end), dlist, 'UniformOutput', 0);
for i = 1 : Ns
    sr = regexprep(crun{i}, '(?i:^run|^ses)([_-])', '');
    sr = strrep(sr, '*', '\w*');
    drun_all{i} = dlist(~cellfun(@isempty, regexp(davr, [sr, '$'], 'match')));
end
drun = unique(vertcat(drun_all{:}));
%%% TO DO - be sure no bug here to find the matching directory...

% Initialize or read run validity
function isval = add_isval(pval)
if ~exist(pval, 'file')
    isval = 1;
    update_valrun(pval, isval);
else
    isval = update_valrun(pval);
end

% Define the general HDR structure holding general hdr info (as output by ft_read_header) 
% + events and headshape:
%- hdr.Fs : sampling frequency in Hz
%- hdr.grad : grad info
%- hdr.event : all events found in raw data as returned by ft_read_event
%- hdr.headshape : headshape as read by ft_read_headshape
function phdr = add_hdr(prawi, pout)
phdr = [pout, filesep, 'hdr_event.mat'];
if exist(phdr, 'file')
    return;
end

% Raw data path
%%% TO DO: see for EEG Biosemi/EGI data format + add the type in hdr_event.type
%%% 'meg' or 'eeg' for further channel selection
praw = filepath_raw(prawi);

% Read general HDR info
hdr_event = ft_read_header(praw);
hdr_event.dataset = praw;

% Add events
hdr_event.event = ft_read_event(praw);

% Add shape if shape file is found (for 4d and neuromag data for now)
pshape = [];
if exist([prawi, filesep, 'hs_file'], 'file')
    pshape = [prawi, filesep, 'hs_file'];
else
    [~, ~, ext] = fileparts(praw);
    if strcmp(ext, '.fif')
        pshape = praw;
    end
end
if ~isempty(pshape)
    hdr_event.headshape = ft_read_headshape(pshape, 'unit', 'mm');
else
    hdr_event.headshape = [];
end

save(phdr, 'hdr_event');

% Make a txt file displaying events and some recording information
cmeg_event_txt(hdr_event, [pout, filesep, 'events_info.txt'], praw);

