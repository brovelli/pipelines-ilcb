function cfold = dirfold(ploc, rmdir)
% Get list of folders at ploc location (only folders and without the parent and
% current folders '..' and '.')
% ploc: path of the folder where to find subfolders
% cfold: list of the folder names found in ploc (cell of string)
% An empty cell is return if no folders were found
% rmdir: optionnal argument: list of folder name to be excluded from cfold
% (additionnally to '..' and '.')
%-CREx-181115
if nargin < 2
    rmdir = [];
else
    if ischar(rmdir)
        rmdir = {rmdir};
    end
end

cfold = {};

dlist = dir(ploc);  
if isempty(dlist)
    return;
end

dnam = {dlist.name}';
% All folders excepted '.' and '..'
isf = [dlist.isdir]' & cellfun(@(x) ~strcmp(x(1),'.'), dnam);
cfold = dnam(isf);
if ~isempty(rmdir)
    cfold = cfold(~ismember(cfold, rmdir));
end
    