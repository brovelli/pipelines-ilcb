function cp_fwd_model_fig(model, pfig)
% Figures of the forward model geometry: sources + conduction volume + MEG sensor
%
% model: forward model as output by cp_fwd_leadfield
% model is a structure with fields 'cortical' and 'subcortical'
%
%-CREx180704
%

if ~isempty(model.subcortical)
    plot_geom(model.subcortical, pfig, 'subcortical')
end

if ~isempty(model.cortical)
    plot_geom(model.cortical, pfig, 'cortical')
end

function plot_geom(modl, pfig, styp)

ins = modl.grid.inside;
dip_pos = modl.grid.pos(ins, :);
lab = modl.atlas.label(ins);

vol = modl.headmodel;
% Define dipole spheres with radius = 0.6 mm

Sballs = dipole_spheres(dip_pos, lab, 0.6, [0.85 0 0]);

Ndip = length(Sballs);

figure
set(gcf, 'visible', 'off', 'color', [0 0 0], 'units', 'centimeters', 'position', [5 5 20 18])
set(gca, 'color', [0 0 0], 'position', [0 0 1 1])
hold on

for n = 1: Ndip
    patch(Sballs(n));
end

ctx_patch = struct('vertices',vol.bnd.pos, 'faces', vol.bnd.tri);

patch(ctx_patch, 'facecolor', 'none',...
                'edgecolor', [0.70    0.98    0.70],...
                'edgealpha', 0.15,...
                'PickableParts', 'none'); 
view(165, 10)
axis equal tight 
set(gca, 'position', [0.02 0.02 0.96 0.8]);
axis off

% dispname function in ft_CREx/atlas toolbox
set(gcf, 'WindowButtonDownFcn', @dispname);

camlight('left','infinite')

% %- Add the famous general title
annotation(gcf,'textbox', 'String', ['Forward model geometry - ', styp],...
    'interpreter','none','FontSize',13,'fontname','AvantGarde','color',[1 1 1],...
    'LineStyle','none','HorizontalAlignment','left',...
    'FitBoxToText','off','Position',[ 0.0172 0.9236 0.9667 0.0661]);
    
ht = title('Click on source to display its atlas region', 'color', [1 1 1], 'fontsize', 14);
export_fig([pfig, filesep, 'forward_model_', styp, '.png'], '-m2');
save_fig([pfig, filesep, 'forward_model_', styp, '.fig']) 

% Add MEG channels
delete(ht)
Sgrad = modl.grid.cfg.grad;
plot_megchan(Sgrad)
view(180, 0)

export_fig([pfig, filesep, 'forward_model_', styp, '_megsens.png'], '-m2');
save_fig([pfig, filesep, 'forward_model_', styp, '_megsens.fig']) 
close 

% Define the sphere patches at dipole locations 
function Sbals = dipole_spheres(dip_pos, lab, radius, col)

Ndip = length(dip_pos);

%----
% Set dipoles spheres
[Xsph, Ysph, Zsph] = sphere(10); 
% Unit sphere with N=42 faces - Xb, Yb and Zb are (N+1) x (N+1) matrices
% using by surf(Xb, Yb, Zb) to draw the sphere
% Set radius of sphere representing dipoles = 6 mm

Xb = Xsph.* radius;
Yb = Ysph.* radius;
Zb = Zsph.* radius;

%----
% Define dipoles patch balls properties

Sbals = struct;
for n = 1: Ndip 

    % Coordinates of the center
    cx = Xb + dip_pos(n, 1);
    cy = Yb + dip_pos(n, 2);
    cz = Zb + dip_pos(n, 3);

    [Sbals(n).faces, Sbals(n).vertices] = surf2patch(cx, cy , cz );   
    Sbals(n).facecolor = col;
    Sbals(n).edgecolor = 'none';
    Sbals(n).facelighting = 'gouraud';
    Sbals(n).displayname = lab{n};
end
