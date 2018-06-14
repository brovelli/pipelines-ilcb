function [pdat, ndat] = find_files(dpath, cext)
if nargin < 2
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

if isempty(pdat)    
    fprintf('Looking for list of files with extension:')
    disp(cext)
    warning('Data files not found in directory\n%s\n', dpath)
end