function Sdb = cp_marsatlas(Sdb)
% Prepare atlas files (surf, vol) 
% - Realign surf and vol according to Mtrans_ref transformation matrix in the
% subject MRI_bv space
% - Realign according to fiducial previously indicated (cp_fwd_singleshell)
%
%   --> marsatlas.mat saved inside db_ft/PROJ/SUBJ/atlas
%       == marsatlas structure with fields "surf" and "vol"
%
%-CREx180530

% Number of data to process
Np = length(Sdb);

for i = 1 : Np
    
    % Check if marsatlas already defined
    if isempty(Sdb(i).anat.atlas)
        Sdb(i) = prepare_atlas(Sdb(i));
    end     
end

function dps = prepare_atlas(dps)
dpa = dps.anat;
%-- Read MarsAtlas surface files from BV/FS pipeline (Left and Right) 
% + textures information
masurf = read_mars_surf(dpa.surf, dpa.tex);
% Number of surf files
Nh = length(masurf);

%-- Read MarsAtlas volume
pvol = dpa.vol{1};
mavol = read_mars_vol(pvol);

%-- Transformation matrices
%* to realign surf and vol in subject MRI space
Meach = loadvar(dpa.trans);
Msu2bv = Meach{1}*Meach{2};
Mvo2bv = inv(Meach{4}*Meach{3});

%* to realign according to fiducial
Mreal = loadvar(dpa.Mreal);    

%-- Realign MarsAtlas

%- Surface
% Total transform
Mt_su = Mreal*Msu2bv;
for k = 1 : Nh
    masurf{k} = ft_transform_geometry(Mt_su, masurf{k});
end

%- Volume
% Total transform
Mt_vo = Mreal*Mvo2bv;  %#ok
% The transform field of mavol atlas will be change according to previous
% transfom + Mt_vol
mavol = ft_transform_geometry(Mt_vo, mavol); 

% Save atlas file
marsatlas = [];
marsatlas.surf = masurf;
marsatlas.vol = mavol;

patlas = [make_dir([dpa, filesep, 'atlas']), filesep, 'marsatlas.mat'];
save(patlas, 'marsatlas');

dps.anat.atlas = patlas;

% Load volume of condution for figure
if isfield(dps.fwd, 'shell') && exist(dps.fwd.shell, 'file')
    shell = loadvar(dps.fwd.shell);
else
    shell = [];
end

% Figure will be saved in coreg directory
pcor = make_dir([dps.dir, filesep, 'coreg']);

% Neeed for raw MEG data file to superimpose sensors
pmeg = dps.meg.continuous.raw{1};

% Check for coregistration of all the things of the universe
coreg_fig(marsatlas, shell, pmeg, pcor);  

% Figure showing superimposition of marsatlas surface + volume, + head model
% with MEG sensor
function coreg_fig(marsatlas, shell, pmeg, pfig)

masurf = marsatlas.surf;
mavol = marsatlas.vol;

% Determine right hemisphere for one of the specific plot
ir = cellfun(@(x) strcmp(x.name, 'surf_R'), masurf);

%-- Visualize

% Head model + surfaces
opt = [];
opt.visible = 'off';
if ~isempty(shell)
    shell.name = 'shellcond';
end

plot_meshes([{shell}; masurf], opt);
view(30, 10)
title('Single-shell conduction volume and MarsAtlas surface')
axis tight
    

% Print
export_fig([pfig, filesep, 'marsatlas_surf_shellcond.png'], '-m2')
save_fig([pfig, filesep, 'marsatlas_surf_shellcond.fig'])
close

% Headmodel, Right MA surface + Left MA volume
plot_vol(mavol, 1:41, 0)

% Superimpose Mesh
opt.newfig = 0;
opt.colors = [0 0.9 0.7; 1 0.4 0.4];
plot_meshes({masurf{ir}, shell}, opt);
view(122, 10)
rotate3d off
title('Click on a MarsAtlas parcel to get its name', 'color', [1 1 1])

% Print
export_fig([pfig, filesep, 'marsatlas_surfR_vol_shellcond.png'], '-m2')
save_fig([pfig, filesep, 'marsatlas_surfR_vol_shellcond.fig'])
close

% Add the MEG sensors
%- Looking for the raw MEG data in pmeg directory
draw = filepath_raw(pmeg);
if isempty(draw)
    return;
end

% Check reg MEG + vol model only
Sgrad = ft_read_sens(draw);

% Everything on the same figure (volumetric MarsAtlas (Right) + the Right
% MarsAtlas surface + the singleshell head model + the MEG sensors
plot_vol(mavol, 1:41, 0);

% surface_R + shell volume
opt.newfig = 0;
[hm, hleg] = plot_meshes({masurf{ir}, shell}, opt);

% MEG sensors
plot_megchan(Sgrad)

view(160, 10)    
title('Cliquer sur une parcelle MarsAtlas pour afficher son nom', 'color', [1 1 1])

%- Print
export_fig([pfig, filesep, 'marsatlas_surfR_vol_shellcond_megsens.png'], '-m2')
save_fig([pfig, filesep, 'marsatlas_surfR_vol_shellcond_megsens.fig'])

% Try to add the headshape
hs = read_headshape(pmeg);

if isempty(hs)
    close
    return;
end

% Headshape (MEG source)
ps = plot3c(hs.pos,'kd');
set(ps, 'markerfacecolor', [0.85 0 0], 'markersize',5, 'PickableParts', 'none');

% Add fiducials
fid = hs.fid;
% Nasion;
isf = cellfun(@(x) strcmpi(x(1), 'n'), fid.label);
pnas = plot3c(fid.pos(isf, :), 'd');
set(pnas, 'markerfacecolor', [0 0.4 0.9], 'markersize', 16,...
    'markeredgecolor', [0.9 0.75 0], 'linewidth', 2);

% LPA
isf = cellfun(@(x) strcmpi(x(1), 'l'), fid.label);
pleft = plot3c(fid.pos(isf, :), 'd');
set(pleft, 'markerfacecolor', [0.9 0.75 0], 'markersize', 16,...
    'markeredgecolor', [0 0.9 0], 'linewidth', 2);

% RPA
isf = cellfun(@(x) strcmpi(x(1), 'r'), fid.label);
pright = plot3c(fid.pos(isf, :), 'd');
set(pright, 'markerfacecolor', [0.85 0.33 0.1], 'markersize', 16,...
    'markeredgecolor', [0 0.75 0.75], 'linewidth', 2);

view(90, 0)

delete(hleg)
put_legend([hm ; pnas ; pleft ; pright], {'surf_R', 'shellcond', 'NAS', 'LPA', 'RPA'});

%- Print
export_fig([pfig, filesep, 'marsatlas_surfR_vol_shellcond_megsens_fid.png'], '-m2')
save_fig([pfig, filesep, 'marsatlas_surfR_vol_shellcond_megsens_fid.fig'])
close

% Add meshes legend
function lg = put_legend(hdl, names) 

% Plot legend 
lg = legend(hdl, names,'location','eastoutside');   

% Set new position
pos = [0.84199 0.60382 0.14519 0.3193]; 
set(lg, 'position', pos,...
    'interpreter', 'none', 'fontsize', 11,...
    'EdgeColor', 'none', 'color', [0.98 0.96 0.96])  

% The color rectangles are of type 'patch' and the
% associated texts of type 'text'.

% Reduced rectangle size
cpa = findobj(lg, 'type', 'patch');
for j = 1 : length(cpa)
    xd = get(cpa(j), 'XData');
    newxd = [xd(1:2); xd(3:4)./2]; 
    set(cpa(j), 'XData', newxd)
end
% Move text closer to the color rectangle
ctx = findobj(lg, 'type', 'text');
for k = 1 : length(ctx)
    pos = get(ctx(k), 'position');
    pos(1) = newxd(3) + 0.02; 
    set(ctx(k),'position', pos);
end 

