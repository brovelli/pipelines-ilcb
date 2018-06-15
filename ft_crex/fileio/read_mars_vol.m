function vol = read_mars_vol(pvol)
% Read marsatlas volumetric file by using ft_read_mri
[pdir, vnam] = fileparts(pvol);

% altas_label & labelinfo from ft_crex toolbox
plabel = which('atlas_label.txt');
pinfo = 'labelinfo.mat';

% Copy atlas_label file in th esame directory as volume file
% Rename atlas_label with volume file name
copyfile(plabel, [pdir, vnam, '.txt']);

labelinfo = loadvar(pinfo, 'labelinfo');

% ft_read_atlas will consider an AAL atlas as nii file is associated with a txt file
vol = read_atlas_ft(pvol, labelinfo);

% Add labelinfo (could be usefull for further parcel search/grouping
vol.labelinfo = labelinfo;

% Read altas - based on ft_read_atlas function (aal case)
% Remove the atlas tissue that are not associated with atlas labels
function atlas = read_atlas_ft(pvol, labelinfo)

lab = labelinfo.label;
idx = labelinfo.index;

% Read MRI with ft_read_mri and check for units (be sure to be as default units
% which are 'mm')
atlas = read_mri(pvol);

atlas.tissue = atlas.anatomy;
atlas = rmfield(atlas, 'anatomy');

atlas.coordsys = 'mni'; % ? don't remember why...
atlas.tissuelabel       = {};
atlas.tissuelabel(idx)  = lab;

% Set tissue to 0 when it is associated with an empty label
tis = atlas.tissue(:);
iemp = find(cellfun(@(x) isempty(x), atlas.tissuelabel)==1);
Ne = length(iemp);
for j = 1 : Ne
    istis = tis==iemp(j);
    if sum(istis) > 0
        tis(istis) = 0;
    end
end
atlas.tissue = reshape(tis, atlas.dim);

% Make the list with labels compact. See ft_read_atlas
[a, ~, k] = unique(atlas.tissue);
atlas.tissue = reshape(k-1, atlas.dim);
atlas.tissuelabel = atlas.tissuelabel(a(a~=0));
