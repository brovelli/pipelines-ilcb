function Sdb = cp_sources_subcort(Sdb, dsources_mm)
% Set subcortical sources from marsatlas volume
%   --> sources.mat saved inside db_ft/PROJ/SUBJ/sources
%       with new field "subcortical" containing sources positions + label info
%
%-CREx180530

Np = length(Sdb);

radius = dsources_mm/2;
vol_sphere = 4/3*pi*radius^3;

for i = 1 : Np
    
    %-- Data paths
    psubj = Sdb(i);
    
    % Check if other sources already process (cf. cortical)
    adone = 0;
    if isfield(psubj, 'sources')
        pso = fileparts(psubj.sources);
        sources = loadvar(psubj.sources);
        if isfield(sources, 'subcortical') && ~isempty(sources, 'subcortical')
            adone = 1;
        end
    else
        pso = make_dir([psubj.dir, filesep, 'sources']);
        sources = [];
    end
    
    if ~adone
        % Load atlas
        atlas = loadvar(psubj.atlas);

        mavol = atlas.vol;

        sources.subcortical = set_subsources(mavol, vol_sphere, pso);

        save([pso, filesep, 'sources.mat'], 'sources');
        Sdb(i).sources = [pso, filesep, 'sources.mat'];
    end

end

function subcor = set_subsources(atlas, vol_sphere, pso)

% Index of subcortical regions in mavol.tissue
isub = find(strcmp(atlas.labelinfo.lobe, 'Subcortical')==1);

% Corresponding labelinfo
labinfo = structfun(@(x) x(isub), atlas.labelinfo, 'UniformOutput', 0);

Ns = length(isub);

% Set sources for each subcortical region

% Get all atlas position in [Np x 3] vector and associated id
[apos, id_tis] = get_atlaspos(atlas);

% Process each subregion
all_pos = cell(Ns, 1);
% Keep all labelinfo 
Slab = init_labinfo(labinfo);

subcor = [];
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
subpos = cell2mat(all_pos);

% Save subcortical sources
subcor.pos = subpos;
subcor.label = Slab;


%-- Figure
% Figure showing subcortical region + sources
plot_vol(atlas, isub);
plot3(subpos(:, 1), subpos(:, 2), subpos(:, 3), 'og')
title('Check for subcortical sources', 'color', [1 1 1])
saveas(gcf, [pso, filesep, 'sources_subcortical.fig'])
export_fig([pso, filesep, 'sources_subcortical.png'], '-m2')
close
    
function Slab = add_labinfo(Slab, subinfo, Nsrc)
[fnam, Nf] = get_names(Slab);

for i = 1 : Nf
    nam = fnam{i};
    subi = subinfo.(nam);
    Slab.(nam) = [Slab.(nam) ; repmat(subi, Nsrc, 1)];
end

function Slab = init_labinfo(labinfo)
[fnam, Nf] = get_names(labinfo);
Slab = [];
for i = 1 : Nf
    Slab.(fnam{i}) = [];
end

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