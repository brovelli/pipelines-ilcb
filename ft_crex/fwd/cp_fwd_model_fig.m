function cp_fwd_model_fig(Spos, vol)


dip_pos = Spos.pos;
lab = Spos.label.label;
Ndip = length(dip_pos);

%----
% Set dipoles spheres
[Xsph, Ysph, Zsph] = sphere(42); 
% Unit sphere with N=42 faces - Xb, Yb and Zb are (N+1) x (N+1) matrices
% using by surf(Xb, Yb, Zb) to draw the sphere
% Set radius of sphere representing dipoles = 6 mm
radius = 0.6;
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
    Sbals(n).facecolor = [0.85 0 0];
    Sbals(n).edgecolor = 'none';
    Sbals(n).facelighting = 'gouraud';
    Sbals(n).displayname = lab{n};
end

figure
set(gcf, 'color', [0 0 0], 'units', 'centimeters', 'position', [5 5 20 18])
set(gca, 'color', [0 0 0], 'position', [0 0 1 1])
hold on
for n = 1: Ndip
    patch(Sbals(n));
end

ctx_patch = struct('vertices',vol.bnd.pos, 'faces', vol.bnd.tri);

patch(ctx_patch, 'facecolor', 'none',...
                'edgecolor', [0.70    0.98    0.70],...
                'edgealpha', 0.15,...
                'PickableParts', 'none'); 
view(127, 6)
axis equal tight 
set(gca, 'position', [0.02 0.02 0.96 0.8]);
axis off

% dispname function in ft_CREx/atlas toolbox
set(gcf, 'WindowButtonDownFcn', @dispname);

camlight('left','infinite')

% %- Add the famous general title
annotation(gcf,'textbox', 'String', 'Forward model geometry',...
    'interpreter','none','FontSize',13,'fontname','AvantGarde','color',[1 1 1],...
    'LineStyle','none','HorizontalAlignment','left',...
    'FitBoxToText','off','Position',[ 0.0172 0.9236 0.9667 0.0661]);
    

