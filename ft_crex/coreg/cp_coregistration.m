function Sdb = cp_coregistration(Sdb)
% Check for coregistration with one of the 2 methods:
%
% - if a headshape file containing fiducial information was found with MEG data,
% distances between fiducial points (LPA, RPA) and the mean right side position of the 
% MarsAtlas surface are compared to determine if a Left/Right flip was done
% during fiducial identification on the MRI image
%
% - if no headshape was found, the coregistration check is done considering the
% normals that are computing on the MarsAtlas surface: they should point outward of
% the surface if the LPA/RPA have been correctly identified, and inward if a
% Left/Right flip has been done during fiducial identification
%
%-CREx181211

% Coregistration check is done when MEG data are ready
imeg = is_meg(Sdb);
if ~any(imeg)
    return
end

Sdbm = Sdb(imeg);
Np = length(Sdbm);

% Initialize waitbar
wb = waitbar(0, 'Check for coregistration...', 'name', 'Coregistration');
wb_custcol(wb, [0 0.6 0.8]);

for i = 1 : Np
    
    waitbar((i-1)/Np, wb, ['Coregistration check: ', Sdbm(i).sinfo]);
    
    dps = Sdbm(i);
    
    if dps.coreg || ~all(dps.meg.run.valid)
        continue;
    end
    
    % Conduction volume
    vol_shell = loadvar(dps.fwd.shell);
    
    % MarsAtlas
    marsatlas = loadvar(dps.anat.atlas);
    
    % FID in headshape
    hdr = loadvar(dps.meg.hdr_event{1});
    
    if any(cellfun(@isempty, {vol_shell, marsatlas, hdr}))
        continue;
    end
    masurf = marsatlas.surf;
    % Right hemisphere
    ir = cellfun(@(x) strcmp(x.name(end), 'R'), masurf);
    surf_r = masurf{ir};
    if ~isempty(hdr.headshape) && isfield(hdr.headshape, 'fid')
        is_flip = check_hs(surf_r, hdr.headshape.fid);
    else
        is_flip = check_nrm(surf_r);
    end
    % Coreg folder
    pcor = make_dir([dps.dir, fsep, 'coreg']);
    % Save coreg.fid_flip indication
    coreg = [];
    coreg.fid_flip = is_flip;
    pmat = [pcor, fsep, 'coreg.mat'];
    save(pmat, 'coreg');
    
    % Flip the marsatlas surface and volume and the singleshell
    if is_flip
        warning(['%s\nLeft/Right fiducial flip has been detected... \n',...
            'Flipping the marsatlas and singleshell volume for source analysis'], dps.sinfo);       
        [marsatlas, vol_shell] = flip_fwd(marsatlas, vol_shell);
        save(dps.anat.atlas, 'marsatlas')
        save(dps.fwd.shell, 'vol_shell');
    end
    
    % Figure will be saved in coreg directory   
    % Check for coregistration of all the things of the universe
    coreg_fig(marsatlas, vol_shell, hdr.grad, hdr.headshape, pcor, dps.sinfo);   
    
    Sdbm(i).coreg = 1;
end
close(wb);
Sdb(imeg) = Sdbm;

% Check if left-right flip when defining fiducial from fid distance form the
% mean right surface position
function isflip = check_hs(surf_r, fid)
mpos = mean(surf_r.pos);
% LPA
isf = cellfun(@(x) strcmpi(x(1), 'l'), fid.label);
lpa = fid.pos(isf, :);
% RPA
isf = cellfun(@(x) strcmpi(x(1), 'r'), fid.label);
rpa = fid.pos(isf, :);

isflip = norm(rpa - mpos) > norm(lpa - mpos);

% Check if left-right flip according to normals inversion on marsatlas surface
% file (a L/R flip will produce inwards normals instead of outwards)
function isflip = check_nrm(surf)
ori = mesh_norm(surf);
% Estimate the volume of the surface mesh
[~, Vsu] = convhull(surf.pos);
% Volume of the surface mesh points shifted from the ori value
% (should be > Vsu if normals are outwards)
[~, Vn] = convhull(surf.pos + ori);
isflip = Vsu > Vn;

% Flip L/R
function [matlas, vol] = flip_fwd(matlas, vol)
% Matrix for flipping (mirrored on y)
Mflip = eye(4, 4);
Mflip(2, 2) = -1;

%- Surface
masurf = matlas.surf;
Nh = length(masurf);
for k = 1 : Nh
    masurf{k} = ft_transform_geometry(Mflip, masurf{k});
end
matlas.surf = masurf;

%- Volume
matlas.vol = ft_transform_geometry(Mflip, matlas.vol); 

%- Shell 
vol = ft_transform_geometry(Mflip, vol);

% Figure showing superimposition of marsatlas surface + volume, + head model
% with MEG sensor
function coreg_fig(marsatlas, shell, grad, hs, pfig, subj)

%-- MEG channels + head conduction volume
cmeg_fwd_checkreg_fig(shell, grad, subj, pfig)

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

% Everything on the same figure (volumetric MarsAtlas (Right) + the Right
% MarsAtlas surface + the singleshell head model + the MEG sensors
plot_vol(mavol, 1:41, 0);

% surface_R + shell volume
opt.newfig = 0;
[hm, hleg] = plot_meshes({masurf{ir}, shell}, opt);

% MEG sensors
plot_megchan(grad)

view(160, 10)    
title('Cliquer sur une parcelle MarsAtlas pour afficher son nom', 'color', [1 1 1])

%- Print
export_fig([pfig, filesep, 'marsatlas_surfR_vol_shellcond_megsens.png'], '-m2')
save_fig([pfig, filesep, 'marsatlas_surfR_vol_shellcond_megsens.fig'])

if isempty(hs)
    close
    return;
end

% Headshape (MEG source)
ps = plot3c(hs.pos,'kd');
set(ps, 'markerfacecolor', [0.85 0 0], 'markersize',5, 'PickableParts', 'none');

if ~isfield(hs, 'fid')
    close
    return;
end
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