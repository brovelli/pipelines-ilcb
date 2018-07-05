function Smeg = cp_db_megpaths(pmeg)
%-- Prepare MEG data paths from the meg directory in FIELDTRIP_DATABASE
% The MEG data directory as prepare by cp_db_import should have this
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
%-CREx-180702

prepare_db_meg(pmeg);

Smeg = [];

pcont = fullfile(pmeg, 'continuous');

% No run directory
if exist([pcont, filesep, 'raw'], 'dir')
    % Add paths to continuous dataset 
    Smeg.continuous.raw = {[pcont, filesep, 'raw']};
    % Add preproc directory if not done yet
    Smeg.continuous.prep = {make_dir([pcont, filesep, 'preproc'])};
    % Add paths to the epoched data
    Smeg.epoched.run = {make_dir(fullfile(pmeg, 'epoched'))};
    Smeg.epoched.concat = Smeg.epoched.run{1};
    return;
end

% Run directories
drun = dir([pcont, filesep, 'run_*']);
if isempty(drun)
    warning('No MEG data found in %s\n', pmeg)
    return
end

Nd = length(drun);
% Initialize data paths
Smeg.continuous.raw = cell(Nd, 1);
Smeg.continuous.prep = cell(Nd, 1);
Smeg.epoched.run = cell(Nd, 1);
for j = 1 : Nd
    rnam = drun(j).name;
    prun = [pcont, filesep, rnam];
    % Add run directories data paths
    Smeg.continuous.raw{j} = [prun, filesep, 'raw'];
    Smeg.continuous.prep{j} = make_dir([prun, filesep, 'preproc']);
    % Add paths to the epoched data
    Smeg.epoched.run{j} = make_dir(fullfile(pmeg, 'epoched', rnam));
    
end

Smeg.epoched.concat = make_dir(fullfile(pmeg, 'epoched', 'run_concat'));


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
        
    % Move all pmeg content in pmeg / continuous / raw directory
    praw = make_dir([pcont, filesep, 'raw']);  
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

