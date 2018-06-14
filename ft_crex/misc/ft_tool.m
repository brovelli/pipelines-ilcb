function statut = ft_tool(toolname, sopt, extft)
% Add toolbox directories (initially to work with fieldtrip processings)
% sopt : option string to indicate if toolbox is to add (sopt = 'add') or to
% remove from matlab's current paths (sopt='rm') [ default: sopt = 'add' ]
% extft flag : = 1 to add toolbox that is external to the "fieltrip/external" directory
%
%-CREx-180411

if nargin < 2 || isempty(sopt) || ~ischar(sopt)
    sopt = 'add';
end
if nargin < 3
    extft = 0;
end


if ~extft 
    % Try to add toolbox with ft_hastoolbox
    % Sometimes cause error => try catch to avoid it
    try 
        statut = ft_hastoolbox(toolname, 1);
    catch
        statut = 0;
    end
    if statut
        % Remove paths
        if strcmp(sopt, 'rm')
            % Get the paths (main + sub-directories)
            ptools = get_added;
            rmpaths(ptools)
        end
        return;
    end
end
    
% Else, try to find this damned toolbox directory !
if strcmp(sopt, 'add')
    ptool = search_tool(toolname, extft);
    if isempty(ptool)
        statut = 0;
        warning(['Impossible to find ', toolname,' toolbox - Check for toolbox directory name'])
    else
        statut = 1;
        fprintf('\nAdding %s toolbox to paths\n', ptool);
        addpath(genpath(ptool))
    end
    return;
end

if strcmp(sopt, 'rm')
    % Find from the Matlab path the previously added toolbox directory
    rmpaths(get_tooldirs(toolname));
end

function rmpaths(ptools)
% Matlab rmpath function removes only one directory with input as path string
if ~isempty(ptools)
    Np = length(ptools);
    for i = 1 : Np
        rmpath(ptools{i})
    end
    fprintf('Removing toolbox and subdirectories \n%s \n', ptools{1});
end

% Find the previously added tool
function ptools = get_tooldirs(toolname)
allpath = strsplit(path, ';');

ismain = cellfun(@(x) strcmpi(x, [fileparts(x), filesep, toolname]), allpath);
if any(ismain)
    pmain = allpath(find(ismain==1, 1, 'first'));
    ptools = find_alldir(pmain, allpath);
else
    ptools = [];
end

% Find the last added toolbox + subdirectories
function ptools = get_added
% Last added toolbox
allpath = strsplit(path, ';');
pmain = allpath{1};
ptools = find_alldir(pmain, allpath);

function ptools = find_alldir(maindir, allpath)
% Find all children directories (inside pmain) + pmain
isins = cellfun(@(x) ~isempty(strfind(x, maindir)), allpath);  %#ok
ptools = allpath(isins);

% Try to find the toolbox directory in main paths (fieldtrip/external,
% fieldtrip, directory holding fieldtrip, and pwd)
function ptool = search_tool(toolname, extft)
if nargin < 2
    extft = 0;
end
pft = which('ft_defaults');
% Find if toolbox is in the same directory where fieldtrip is 
ftdir = fileparts(pft);
extdir = [ftdir, filesep, 'external'];

% Matlab toolbox directory (where filedtrip is installed)
mtool = fileparts(ftdir);

% Matlab toolbox directory
mmain = toolboxdir('');
if extft
    spath = {mmain ; mtool ; pwd};
else
    spath = {extdir ; ftdir; mtool ; mmain ;  pwd};
end
Np = length(spath);
ptool = [];
for i = 1 : Np
    sdir = spath{i};
    allp = strsplit(genpath(sdir), ';');
    isdir = strcmpi(allp, [sdir, filesep, toolname]);
    if any(isdir)
        ptool = allp{find(isdir==1, 1, 'first')};
        break;
    end
end
