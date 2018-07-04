function fpath = make_dir(dirpath, new)
% Make a directory to store figures / results
% dirpath : path of the directory to create
% new : flag indicating if a prefix is to add to the new directory name (an
% incremential number depending on the last maximum number found in any pre-existing
% directories (nmax). Final name of the new directory will be : [dirpath, '_', nmax + 1]
%
% CREx-20170613

if nargin < 2 || isempty(new)
    new = 0;
end

% Directory already exists
if ~new && exist(dirpath, 'dir')
    % Return the full path
    fpath = fullpath(dirpath);
    return
end

% Remove double filesep 
dirpath = regexprep(dirpath, '\\\', '\\');

% Remove first filesep (path at working directory)
if strcmp(dirpath(1), filesep)
    dirpath = dirpath(2:end);
end

% Remove last filesep
if strcmp(dirpath(end), filesep)
    dirpath = dirpath(1:end-1);
end

% Check if other directories already here
if new==1
    % Directory name will hold a suffix with directory number
    % to ensure not to use previous directory
    dirpath = [dirpath, '_'];
    ocheck = dir([dirpath, '*']);
else
    ocheck = dir(dirpath);
end

% Create the new directory     
if isempty(ocheck)
    if new
        dirpath = [dirpath, '1'];
    end
    mkdir(dirpath);
else
    % Be sure the prefix match with the directory name and it's following only
    % by the suffix_* number
    if new
        [~, ininam] = fileparts(dirpath);
        cnam = {ocheck(:).name};
        cnum = regexp(cnam, ['(?:', ininam(1:end-1), ')_(\d+)$'], 'tokens');
        inum = cellfun(@(x) ~isempty(x), cnum);
        if any(inum)
            dirnum = cellfun(@(x) str2double(x{1}{1}), cnum(inum), 'uniformoutput', 1);
            nmax = max(dirnum);
            dirpath = [dirpath, num2str(nmax + 1)];
        else
            dirpath = [dirpath, '1'];
        end
        
        mkdir(dirpath);
    end
end
fpath = fullpath(dirpath);

% Return the full path
function dirpath = fullpath(dirpath)
if isempty(strfind(dirpath, filesep))  %#ok
    dirpath = [pwd, filesep, dirpath];
end