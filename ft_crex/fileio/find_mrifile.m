function [pmri, nmri] = find_mrifile(dirpath, imgext)
if nargin < 2
    imgext = {'mri','nii'};
end
if ~iscell(imgext)
    imgext = {imgext};
end
Nf = length(imgext);

for i = 1:Nf
    sext = imgext{i};
    [pmri, nmri] = dirlate(dirpath, ['*.', sext]);
    if ~isempty(pmri)
        break;
    end
end

if isempty(pmri)
    fprintf('MRI file not found in directory\n%s\n', dirpath)
    fprintf('Check for file formats testing in this code: find_mrifile.m\n\n')
end