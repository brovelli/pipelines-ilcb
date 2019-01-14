function hfig = dynmesh_plot(Smesh, pow, opt)
% Figures showing the mesh cortical activity with 4 subplots showing 4 differents views
% according to opt.view

%-- Check for options

% Default width for each subplot
dw = [0.45 0.45];
% Default subplot positions ('pos'), associated views ('views')
dopt = struct('pos', [0.05 0.46 dw ; 0.54 0.46 dw ; 0.05 0 dw ; 0.54 0 dw],...
            'views', [-90 90 ; 90 -90 ; -180 0 ; 0 0],...
            'colmap', [],...
            'clim', [],...
            'title', '',...
            'timer', '',...
            'condition', '');

opt = check_opt(opt, dopt);

if isempty(opt.clim)
    opt.clim = [min(pow(:)) max(pow(:))];
end
if isfield(Smesh, 'norm')
    ori = Smesh.norm;
    Smesh = rmfield(Smesh, 'norm');
else
    ori = [];
end

% Prepare mesh for patch function (remove additionnal fieldnames that will make patch error)
kf = {'faces', 'vertices'};
ik = setdiff(fieldnames(Smesh), kf);
if ~isempty(ik)
    Smesh = rmfield(Smesh, ik);
end
% Make dynamic mesh comparison figure   
hfig = figure('visible', 'off', 'units', 'centimeters', 'position', [2 1 22 22], 'color', [0 0 0]);

pos = opt.pos;
axa = zeros(4, 1);

%- 4 meshes views (from the top, bottom, left and right)
for j = 1 : 4    
    ax = axes(hfig, 'position', pos(j, :));
    hm = patch(Smesh, 'EdgeColor', 'none');
    set(hm, 'FaceColor', [0.38 0.38 0.38], 'FaceAlpha', 0.3)
    
    hs = patch(Smesh, 'EdgeColor', 'none',...
        'FaceVertexCData', pow,...
        'FaceColor', 'interp',...
        'FaceAlpha', 0.95);
    if ~isempty(ori)
        set(hs, 'VertexNormals', ori);
    end
    lighting gouraud
    set(hs, 'AmbientStrength', 0.5, 'DiffuseStrength', 0.8,...
        'SpecularStrength', 0.1, 'SpecularExponent', 8, 'SpecularColorReflectance', 0.1);
    set(ax, 'view', opt.views(j, :))  
    caxis(ax, opt.clim)
    camlight('right', 'infinite')
    axa(j) = ax;
end
axis(axa, 'tight', 'equal', 'vis3d', 'off')

figure(hfig);

%- Change colormap
if ~isempty(opt.colmap)
    colormap(opt.colmap)
    pause(0.002)
end

%- Colorbar
cb = colorbar('color', [1 1 1], 'fontsize', 10, 'FontWeight', 'bold');
set(cb, 'position', [0.515 0.35 0.015 0.21])
set(cb.Label, 'string', 'z-score', 'fontsize', 12, 'color', [1 1 1],...
    'Units', 'normalized',...
    'Position', [2.7 0.5 0], 'FontWeight', 'normal')
%- Timer
annotation(gcf, 'textbox', [0.31 0.84 0.42 0.08], 'color', [1 1 1],...
    'fontweight', 'bold', 'EdgeColor', 'none',...
    'fontsize', 18, 'String', opt.timer, 'FontName', 'Monospaced',...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')

%- General title
annotation(gcf, 'textbox', [0.010114 0.92629 0.97896 0.068003], 'color', [1 1 1],...
    'fontweight', 'bold', 'EdgeColor', 'none',...
    'fontsize', 16, 'String', opt.title, 'interpreter', 'none',...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')

%- Condition information
if ~isempty(opt.condition)
    annotation(gcf, 'textbox', [0.22142 0.01083 0.59446 0.054151], 'color', [1 1 1],...
        'FontAngle', 'italic', 'EdgeColor', 'none',...
        'fontsize', 14, 'String', ['Condition: ', opt.condition], 'interpreter', 'none',...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
end