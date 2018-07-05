function plist = dirpaths(pdir)
% List of file/folders full paths inside pdir directory (without parents directories)
dd = dir(pdir);
if isempty(dd)
    plist = [];
    return
end
% Get the full path
pdir = dd(1).folder;
Nd = length(dd);
plist = cell(Nd, 1);
for i = 1 : Nd
    dnam = dd(i).name;
    if ~strcmp(dnam, '.') && ~strcmp(dnam, '..')
        plist{i} = [pdir, filesep, dnam];
    end
end
plist = plist(~cellfun( @isempty, plist));