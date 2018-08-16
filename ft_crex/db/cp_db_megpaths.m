function Smeg = cp_db_megpaths(pmeg)
%-- Prepare MEG data paths from the meg directory in FIELDTRIP_DATABASE
% The MEG data directories as prepared by cp_db_import should have this
% structure:
% meg / 
%   --- continuous 
%           --- (run_*)
%               --- raw
%               --- preproc
%   --- epoched
%           --- (run_* and run_concat)
%
% If the MEG recording isn't made from several sessions (runs) per
% subject, data (raw and mat) should be directly found in 'continuous/raw', 
% 'continuous/preproc' and 'epoched' directories
%
% The hdr_event structure is created inside each raw data directory in hdr_event.mat file.
% This structure holds general information from the raw data as returned by ft_read_header, 
% the events read by ft_read_event and the headshape if a headshape file was found.
%
%-CREx-180702

prepare_db_meg(pmeg);

Smeg = [];

pcont = fullfile(pmeg, 'continuous');

% Run directories
drun = dir([pcont, filesep, 'run_*']);
if isempty(drun)
    warning('No MEG data found in %s\n', pmeg)
    return
end

Nd = length(drun);
% Initialize data paths
Smeg.rundir = {drun(:).name}';
Smeg.Nrun = Nd;
Smeg.continuous.raw = cell(Nd, 1);
Smeg.continuous.prep = cell(Nd, 1);
Smeg.epoched.run_dir = cell(Nd, 1);
Smeg.clean_dir = cell(Nd, 1);
Smeg.clean_mat = cell(Nd, 1);
for j = 1 : Nd
    rnam = drun(j).name;
    prun = [pcont, filesep, rnam];
    praw = [prun, filesep, 'raw'];
    % Add hdr info
    add_hdr(praw);
    % Add run directories data paths
    Smeg.continuous.raw{j} = praw;
    Smeg.continuous.prep{j} = make_dir([prun, filesep, 'preproc']);
    % Add paths to the epoched data
    Smeg.epoched.run_dir{j} = make_dir(fullfile(pmeg, 'epoched', rnam));    
end

% To see later if necessary (depending on the data processing pipeline)
% Smeg.epoched.concat = make_dir(fullfile(pmeg, 'epoched', 'run_concat'));

% Check if continuous folder exists
% If not, try to find raw data in pmeg directory or in subdirectories that will
% be considered as run folders
function prepare_db_meg(pmeg)

pcont = fullfile(pmeg, 'continuous');

% If 'continuous' folder exists inside pmeg directory, the MEG data importation
% is considered to have been already done by cp_db_import
if exist(pcont, 'dir')
    return
end

% Get initial list of files and folders inside pmeg
pdat = dirpaths(pmeg);
Nd = length(pdat);

% Create continuous directory in pmeg folder
pcont = make_dir(pcont);

% If meg data have been already copied inside pmeg directory but not by
% cp_db_import
% Raw MEG data with no run
if ~isempty(filepath_raw(pmeg))
        
    % Move all pmeg content in pmeg / continuous / run_1 /raw directory
    praw = make_dir(fullfile(pcont, 'run_1', 'raw'));  
    for i = 1 : Nd
        movefile(pdat{i}, praw);
    end
    return
end

% If subdirectories are found containing raw MEG data, these
% directories are assumed to be raw run directories
% --> move it to continuous / run_* / raw directory
isrun = zeros(Nd, 1);
for k = 1 : Nd
    psub = pdat{k};
    % Raw data are found inside psub subdirectory of DB_FT/meg folder
    if isfolder(psub) && ~isempty(filepath_raw(psub))
        isrun(k) = 1;
    end
end
isrun = logical(isrun);
% Only one run
if sum(isrun)==1
    prun = pdat{isrun};
    
    movefile(prun, [pmeg, filesep, 'raw']);
    movefile([pmeg, filesep, 'raw'], [pcont, filesep]);
end
% Manage multi-run recording
if sum(isrun) > 1
    prun = pdat(isrun);
    Nr = length(prun);
    for k = 1 : Nr
        psub = prun{k};
        % Change folder name
        % Folder name
        [~, dnam] = fileparts(psub);
        % Rename it to run_*
        if length(dnam) < 3 || ~strcmp(dnam(1:3), 'run')
            rnam = ['run_', dnam];
            movefile(psub, [pmeg, filesep, rnam]);
            psub = [pmeg, filesep, rnam];
        end
        % Get all files inside the run directory
        pdin = dirpaths(psub);
        % Move all inside raw subdirectory
        praw = make_dir([psub, filesep, 'raw']);
        Nd = length(pdin);
        for j = 1 : Nd
            movefile(pdin{j}, [praw, filesep]);
        end
        % Move run directory inside continuous directory
        movefile(psub, [pcont, filesep]);
    end
end

% Define the general HDR structure holding general hdr info (as output by ft_read_header) 
% + events and headshape:
%- hdr.Fs : sampling frequency in Hz
%- hdr.grad : grad info
%- hdr.event : all events found in raw data as returned by ft_read_event
%- hdr.headshape : headshape as read by ft_read_headshape

function add_hdr(pdir)
phdr = [pdir, filesep, 'hdr_event.mat'];
if exist(phdr, 'file')
    return;
end

% Raw data path
praw = filepath_raw(pdir);

% Read general HDR info
hdr_event = ft_read_header(praw);
hdr_event.dataset = praw;

% Add events
hdr_event.event = ft_read_event(praw);

% Add shape if shape file is found (for 4d and neuromag data for now)
pshape = [];
if exist([pdir, filesep, 'hs_file'], 'file')
    pshape = [pdir, filesep, 'hs_file'];
else
    [~, ~, ext] = fileparts(praw);
    if strcmp(ext, '.fif')
        pshape = praw;
    end
end
if ~isempty(pshape)
    hdr_event.headshape = ft_read_headshape(pshape);
end

save(phdr, 'hdr_event');
