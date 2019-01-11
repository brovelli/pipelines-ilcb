function Sdb = cp_fwd_sources(Sdb, dsources_mm)
% Set cortical and subcortical sources 
%   --> sources.mat saved inside db_ft/PROJ/SUBJ/sources
%       with new field "subcortical" and "cortical" containing both sources positions
%       + label info + orientation for cortical ("mom" field)
%
%- Subcortical sources are defined from the marsatlas volume that was
%   previously defined with cp_marsatlas (marsatlas.mat <-- marsatlas.vol data)
%
%- Cortical sources are defined from the cortical surface that was previously
%   prepared by cp_marsatlas (marsatlas.mat <-- marsatlas.surf data)
%   The normals are computed by patchnormals function (see
%   ft_crex/external/patch_normal).
%
%   TO DO : decimate full resolution surface file if not done by Brainvisa
%   (according to input 'npts_decim' for example)
%
%-CREx180530

Np = length(Sdb);

radius = dsources_mm/2;
vol_sphere = 4/3*pi*radius^3;

% Initialize waitbar
wb = waitbar(0, 'Source preparation...', 'name', 'Forward model');
wb_custcol(wb, [0 0.6 0.8]);

for i = 1 : Np
    waitbar((i-1)/Np, wb, ['Sources: ', Sdb(i).sinfo]);
    %-- Data paths
    psubj = Sdb(i);
    
    % Check if other sources already process (cf. cortical or subcortical fields)
    % Be sure coregistration check was done before source preparation
    if exist(psubj.fwd.sources, 'file') || ~psubj.coreg
        continue;
    end
    
    pso = make_dir([psubj.dir, filesep, 'fwd']);
         
    % Load atlas
    atlas = loadvar(psubj.anat.atlas);
    
    sources = [];
    sources.subcortical = set_subsources(atlas.vol, vol_sphere, pso);
    sources.cortical = set_cortsources(atlas.surf, pso);

    save([pso, filesep, 'sources.mat'], 'sources');

    % Save the sources
    Sdb(i).fwd.sources = [pso, filesep, 'sources.mat'];

    % Figure showing everything
    % Coreg fig showing sources + shell conduction + sensor are saved in coreg directory 
    % while the one showing only the subcortical + cortical sources is saved
    % in fwd directory
    pcor = make_dir([psubj.dir, filesep, 'coreg']);
    fwd_coreg_fig(sources, psubj.fwd.shell, psubj.meg, pso, pcor);
end
close(wb);

% Only cortical sources with an associated MarsAtlas label should be scanning
% while beamforming process -- this is done by setting the 'inside' vector with
% dimension [Npos x 3] with value == false when a source position is not
% associated with an atlas label (cf. for voxel with tex index ==0 or 255 (outside/white
% matter))
function grid = set_cortsources(Csurf, pfig)

if isempty(Csurf)
    warning('MarsAtlas surface data are missing: cortical sources will not be defined')
    grid = [];
    return
end
% Load labelinfo 
pinfo = fullfile(ptool, 'atlas', 'labelinfo.mat');
labelinfo = loadvar(pinfo, 'labelinfo');
ilab = labelinfo.index;

Nh = length(Csurf);

% Compute surface normals for each surface (each hemisphere) and 
% the mesh faces for merge surface
for i = 1 : Nh
    surfh = Csurf{i};
    
    % Add normals
    surfh.mom = mesh_norm(surfh);

    Csurf{i} = surfh;
    % Associated labels - remove 0 and 255 = those who are not in
    % labelinfo.index
    idx = surfh.tex;
    ikeep = ismember(idx, ilab);
    
    surfh.inside = ikeep;
    
    surfh.mtri = surfh.tri;
    % Increment face index for the second hemisphere 
    if i > 1
        ptri = Csurf{i-1}.mtri;
        surfh.mtri = surfh.mtri + max(ptri(:));
    end        
    Csurf{i} = surfh; 
end

% Merge the 2 hemispheres
msurf = merge_surf(Csurf);

alltex = msurf.tex;
% Keep surface name for figure
cnam = msurf.name;

msurf.tri = msurf.mtri;
grid = rmfield(msurf, {'mtri', 'name'});
grid.unit = 'mm';

% Associate the labelinfo with the same number of
grid.label = cortical_label(labelinfo, alltex);

% Figure showing cortical sources on surface
cortical_fig(Csurf, cnam, grid, pfig)


function cortical_fig(Csurf, cnam, grid, pfig)

% All source positions
pos = grid.pos(grid.inside, :);

% Figure showing meshes + normals
col = color_group(4);

opt = [];
opt.names = cnam;
opt.colors = col([2 3], :);
opt.dispnorm = 1;
opt.visible = 'off';
plot_meshes(Csurf, opt);

hq = findobj(gca, 'type', 'quiver');
set(hq, 'AutoScaleFactor', 1.1, 'color', [0.87 0.49 0]);
hold on

% Add sources
hp = plot3(pos(:, 1), pos(:, 2), pos(:, 3), 'o');
set(hp, 'markersize', 2, 'markeredgecolor', [0 0.45 0.74], 'LineWidth', 0.3)
view(90, 0)
title('Cortical sources location and orientation')

% Print
export_fig([pfig, filesep, 'sources_cortical_ori.png'], '-m2')
save_fig([pfig, filesep, 'sources_cortical_ori.fig'])
close

% With labels
plot_cort_label(grid)

% Print
export_fig([pfig, filesep, 'sources_cortical_label.png'], '-m2')
save_fig([pfig, filesep, 'sources_cortical_label.fig'])
close

% Get the corresponding labelinfo to the texture vector alltex
function Slab = cortical_label(labelinfo, alltex)

ilab = labelinfo.index;
Slab = init_labinfo(labelinfo, length(alltex));

fnames = fieldnames(Slab);
Nf = length(fnames);

uid = unique(alltex);
Nd = length(uid);
for j = 1 : Nd
    ind = uid(j);
    
    % Corresponding label 
    if any(ilab==ind)
        subinfo = structfun(@(x) x(ilab==ind), labelinfo, 'UniformOutput', 0);
        aind = find(alltex==ind);
        for k = 1 : Nf
            nam = fnames{k};
            Slab.(nam)(aind, :) = repmat(subinfo.(nam), length(aind), 1);
        end
    end
end
% Merge Left and Right meshes
function msurf = merge_surf(Csurf)

Nh = length(Csurf);
fsu = fieldnames(Csurf{1});
Nf = length(fsu);
msurf = [];
for i = 1 : Nf
    cdat = cell(Nh, 1);
    fnam = fsu{i};
    for j = 1 : Nh
        cdat{j} = Csurf{j}.(fnam);
    end
    if ~ischar(cdat{j})
        msurf.(fnam) = cell2mat(cdat);
    else
        msurf.(fnam) = cdat;
    end
end
    
% Define the subcortical sources by using kmeans (euclidean distance) method
function subcor = set_subsources(atlas, vol_sphere, pfig)

if isempty(atlas)
    warning('MarsAtlas volume data is missing: subcortical sources will not be defined')
    subcor = [];
    return
end
% Index of subcortical regions in mavol.tissue
isub = find(strcmp(atlas.labelinfo.lobe, 'Subcortical')==1);

% Number of subcortical parcels
Ns = length(isub);

% Corresponding labelinfo
labinfo = structfun(@(x) x(isub), atlas.labelinfo, 'UniformOutput', 0);

%-- Set sources for each subcortical region

% Get all atlas position in [Np x 3] vector and associated id
[apos, id_tis] = get_atlaspos(atlas);

% Process each subregion
all_pos = cell(Ns, 1);

% Keep all labelinfo in Slab 
% Initialize Slab
Slab = init_labinfo(labinfo);

%-- Set sources by kmean method at each atlas region
for j = 1 : Ns
    idx = isub(j);

    pos = apos(id_tis==idx, :);

    % Approximate volume in mm^3
    Vsu = length(pos);

    % Approximate number of sources to position in the subvolume
    Nsrc = round(Vsu/vol_sphere);

    % Position sources using kmeans with default distance method (euclidian)
    [~, spos] = kmeans(pos, Nsrc, 'distance', 'sqeuclidean', 'start', 'uniform');

    % Store all pos
    all_pos{j} = spos;

    % Label info
    subinfo = structfun(@(x) x(j), labinfo, 'UniformOutput', 0);
    Slab = add_labinfo(Slab, subinfo, Nsrc);        
end

% Get all source positions in a [ Npos x 3 ] matrix
subpos = cell2mat(all_pos);

% Set subcortical sources
subcor = [];
subcor.pos = subpos;
subcor.label = Slab;
subcor.unit = 'mm';

%-- Figure
% Figure showing subcortical region + sources
plot_vol(atlas, isub);
plot3c(subpos, 'og')
title('Check for subcortical sources', 'color', [1 1 1])
export_fig([pfig, filesep, 'sources_subcortical.png'], '-m2')
save_fig([pfig, filesep, 'sources_subcortical.fig'])
close
 
% Add label info for each source from subcortical region
function Slab = add_labinfo(Slab, subinfo, Nsrc)
[fnam, Nf] = get_names(Slab);

for i = 1 : Nf
    nam = fnam{i};
    subi = subinfo.(nam);
    Slab.(nam) = [Slab.(nam) ; repmat(subi, Nsrc, 1)];
end

% Initialize the labelinfo Slab structure
function Slab = init_labinfo(labinfo, Ni)
if nargin < 2
    Ni = [];
end
[fnam, Nf] = get_names(labinfo);
Slab = [];
for i = 1 : Nf
    nam = fnam{i};
    if isempty(Ni)
        Slab.(nam) = [];
    else
        dati = labinfo.(nam);
        ncol = length(dati(1, :));
        if iscell(dati)
            Slab.(nam) = repmat({''}, Ni, ncol); 
        else
            Slab.(nam) = zeros(Ni, ncol);
        end
    end
end

% Get atlas position as a Np x 3 matrix from a 3D volume
function [pos, id_tis] = get_atlaspos(atlas)

%---- Define atlas coordinates as a N*3 matrix (N = prod(dim) = 902629) 
dim = atlas.dim;
[X, Y, Z]  = ndgrid(1:dim(1), 1:dim(2), 1:dim(3));
Apos   = [X(:) Y(:) Z(:)];

%---- Apply transform to put atlas position in head coordinates
% Left and Right become like fieldtrip and BS convention
Mtr = atlas.transform; 
pos = Apos;
pos(:,4) = 1;
pos = pos * Mtr';
pos = pos(: , 1:3);

%---- Associate tissue identification number
id_tis = atlas.tissue(:);

% Figure showing cortical sources with label
function plot_cort_label(grid, fnew)
if nargin < 2
    fnew = 1;
end
if fnew
    figure
    set(gcf, 'visible', 'off', 'units', 'centimeters', 'Position', [2 2 21 23], 'color', [0 0 0])
    set(gca, 'position', [0.079 0.104 0.788 0.733], 'color', [0 0 0])
end

% Surface
patch('faces', grid.tri, 'vertices', grid.pos, 'edgecolor','none',...
    'facecolor', [0.87 1 0.87], 'facealpha', 0.15, 'facelighting','gouraud', 'pickableparts', 'none'); 
    
axis off

hold on

% Add the dipole source + orientation + labels
ins = grid.inside;

gpos = grid.pos(ins, :);
glob = grid.label.lobe(ins);
glab = grid.label.label(ins);

plot_dip_label(gpos, glab, glob);

% Add norms
gori = grid.mom(ins, :);
hq = quiver3c(gpos, gori);
set(hq, 'AutoScaleFactor', 0.8, 'linewidth', 0.5, 'color', [0.45 0.45 0.35], 'PickableParts', 'none');

view(90, 90)

axis tight equal off;
rotate3d off
set(gcf, 'WindowButtonDownFcn', @dispname);

title('Click on a source to get the MarsAtlas parcel name', 'color', [1 1 1], 'fontsize', 16)

% Plot the sources (dipole) with associated label name
function plot_dip_label(gpos, glab, glob)
ulab = unique(glab);
Nlab = length(ulab);

% One color per label
col = color_group(Nlab);
for j = 1 : Nlab 
    slab = ulab{j};
    ipos = strcmp(glab, slab);
    % Associated lobe
    slob = glob{find(ipos==1, 1, 'first')};
    
    % Associated sources
    ppos = gpos(ipos, :);
    
    plot3c(ppos,'o','markersize',3,...
        'markerfacecolor',col(j,:),'markeredgecolor','none',...
        'displayname', [slob, ' - ', slab])
end

%-- Coreg figure with all sources + volume of conduction + brain surface + sensors
function fwd_coreg_fig(sources, pshell, pmeg, pso, pcor)

if ~isempty(sources.cortical)
    % Plot cortical
    plot_cort_label(sources.cortical)
end

if ~isempty(sources.subcortical)
    % Plot subcortical
    subso = sources.subcortical;
    plot_dip_label(subso.pos, subso.label.label, subso.label.lobe)
    ppos = subso.pos;
    % Highlight subcortical
    pp = plot3c(ppos,'o','markersize',5,'linewidth', 0.8,...
        'markerfacecolor','none','markeredgecolor',[0.9 0.9 0.9]);
    set(pp, 'PickableParts', 'none');
end

view(90, 0)
export_fig([pso, filesep, 'sources_all_label.png'], '-m2')
save_fig([pso, filesep, 'sources_all_label.fig'])

% Add shell
if ~isempty(pshell)
    vol = loadvar(pshell);
     
    opt = [];
    opt.newfig = 0;
	opt.names = 'Singleshell';
    plot_meshes(vol.bnd, opt);
    
    % Print
    export_fig([pcor, filesep, 'fwd_vol_sources_label.png'], '-m2')
	save_fig([pcor, filesep, 'fwd_vol_sources_label.fig'])
end

% Add MEG sensors if pmeg is valid
if isempty(pmeg)
    return;
end
hdr = loadvar(pmeg.hdr_event{1});
if ~isempty(hdr)
    % Check reg MEG + vol model only
    plot_megchan(hdr.grad)
    delete(findobj(gcf, 'type', 'legend'))

    % Print
    export_fig([pcor, filesep, 'fwd_vol_sources_label_sens.png'], '-m2')
    save_fig([pcor, filesep, 'fwd_vol_sources_label_sens.fig'])
end

close 
