function [pdat, ndat] = find_files(dpath, cext, wmsg)
% Return the list of the files with specific format found in a directory
% 
%--- Input
% - dpath : path of the directory
%
% - cext : list of file extensions to look for in dpath directory 
%   The data files path are returned when a matching is found with one of the
%   cext format - cext format being tested one after the other until a matching
%   is found
%   [ default is set for anatomical files:  {'gii', 'nii', 'gz', 'mri'}]
%
% - wmsg : flag to display a warning message if data files not found 
%   [ default: 1 - msg is displayed ]
%
%--- Output 
% - pdat : list of data file paths matching with one of the format in cext
% - ndat : data file names (without the extension)
%
%-CREx180530

% Warning message flag
if nargin < 3
    wmsg = 1;
end

% List of file extensions to find
if nargin < 2 || isempty(cext)
    cext = {'gii', 'nii', 'gz', 'mri'};
end

if ~iscell(cext)
    cext = {cext};
end
Nf = length(cext);

pdat = [];
ndat = [];
for i = 1 : Nf
    sext = cext{i};
    sdd = dir([dpath, filesep, '*.', sext]);

    if ~isempty(sdd)    
        ndat = {sdd(:).name}';
        pdat = cellfun(@(x) [dpath, filesep, x], ndat, 'UniformOutput', 0);
        break;
    end
end

if wmsg && isempty(pdat)    
    fprintf('Looking for list of files with extension:')
    disp(cext)
    warning('Data files not found in directory\n%s\n', dpath)
end