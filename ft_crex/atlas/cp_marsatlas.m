function Sdb = cp_marsatlas(Sdb)
% Prepare atlas files (surf, vol) 
% - Realign surf and vol according to Mtrans_ref transformation matrix in the
% subject MRI_bv space
% - Realign according to fiducial previously indicated (cp_fwd_singleshell)
%
%   --> marsatlas.mat saved inside db_ft/PROJ/SUBJ/atlas
%
%-CREx180530

% Number of data to process
Np = length(Sdb);

for i = 1 : Np
    
    % Check if marsatlas already defined
    if isempty(Sdb(i).atlas)
        Sdb(i) = prepare_atlas(Sdb(i));
    end     
end

function dps = prepare_atlas(dps)

%-- Read MarsAtlas surface files from BV/FS pipeline (Left and Right) 
% + textures information
masurf = read_mars_surf(dps.surf, dps.tex);
% Number of surf files
Nh = length(masurf);

%-- Read MarsAtlas volume
pvol = dps.vol{1};
mavol = read_mars_vol(pvol);

%-- Transformation matrices
%* to realign surf and vol in subject MRI space
Meach = loadvar(dps.trans);
Msu2bv = Meach{1}*Meach{2};
Mvo2bv = inv(Meach{4}*Meach{3});

%* to realign according to fiducial
Mreal = loadvar(dps.Mreal);    

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

patlas = [make_dir([dps.dir, filesep, 'atlas']), filesep, 'marsatlas.mat'];
save(patlas, 'marsatlas');

dps.atlas = patlas;

% Load volume of condution for figure
shell = loadvar(dps.shell);

% Figure will be saved in coreg directory
pcor = make_dir([dps.dir, filesep, 'coreg']);

% Neeed for raw MEG data file to superimpose sensors
pmeg = dps.meg;

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
plot_meshes([{shell}; masurf], {'shellcond'; masurf{1}.name;  masurf{2}.name})
view(30, 10)
title('Single-shell conduction volume and MarsAtlas surface')
axis tight
saveas(gcf, [pfig, filesep, 'marsatlas_surf_shellcond.fig'])
export_fig([pfig, filesep, 'marsatlas_surf_shellcond.png'], '-m2')
close

% Headmodel, Right MA surface + Left MA volume
plot_vol(mavol, 1:41)
plot_meshes({masurf{ir}, shell}, {masurf{ir}.name, 'shellcond'}, [0 0.9 0.7; 1 0.4 0.4], 0)
view(122, 10)
rotate3d off
title('Click on a MarsAtlas parcel to get its name', 'color', [1 1 1])
saveas(gcf, [pfig, filesep, 'marsatlas_surfR_vol_shellcond.fig'])
export_fig([pfig, filesep, 'marsatlas_surfR_vol_shellcond.png'], '-m2')
close

% Add the MEG sensors
%- Looking for the raw MEG data in pmeg directory
draw = filepath_raw(pmeg);
if ~isempty(draw)
    % Check reg MEG + vol model only
    Sgrad = ft_read_sens(draw);

    % Everything on the same figure (volumetric MarsAtlas (Right) + the Right
    % MarsAtlas surface + the singleshell head model + the MEG sensors
    plot_vol(mavol, 1:41)
    plot_meshes({masurf{ir}, shell}, {masurf{ir}.name, 'shellcond'}, [0 0.9 0.7; 1 0.4 0.4], 0)
    plot_megchan(Sgrad)
    view(160, 10)    
    title('Cliquer sur une parcelle MarsAtlas pour afficher son nom', 'color', [1 1 1])
    saveas(gcf, [pfig, filesep, 'marsatlas_surfR_vol_shellcond_megsens.fig'])
    export_fig([pfig, filesep, 'marsatlas_surfR_vol_shellcond_megsens.png'], '-m2')
    close
end