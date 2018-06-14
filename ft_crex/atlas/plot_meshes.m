function plot_meshes(meshes, names, colors, newfig)
% Plot one or several meshes on the same figure
% meshes : cell of meshes or unique mesh data [REQUIRED]
% names : associated names (cell) [default : mesh_1, mesh_2...]
% colors : associated colors as a [Nmesh x 3] matrix [default : color defined by
% color_group function]
%
% Use ft_CREx color_group function 
%
%-CREx-20180206

if ~iscell(meshes)
    meshes = {meshes};
end
Nm = length(meshes);
if nargin < 4
    newfig = 1;
end

if nargin < 3 || isempty(colors)
    colors = color_group(Nm);
end

if nargin < 2 || isempty(names)
    names = default_names(Nm);       
end
if ischar(names)
    names = {names};
end

if newfig
    figure
    set(gcf, 'units', 'centimeters', 'Position', [2 2 21 23])
    set(gca, 'position', [0.079 0.104 0.674 0.7335])
end
hold on
hp = zeros(Nm, 1);
for i = 1 : Nm
    mesh = prep_patch(meshes{i});
    hp(i) = patch(mesh, 'edgecolor', colors(i, :), 'edgealpha', 0.2, 'facecolor', [0.95 0.90 0.80],...
        'facealpha', 0.15, 'facelighting','gouraud', 'pickableparts', 'none'); 
end
put_legend(hp, names);

xlabel('x', 'fontsize', 14)
ylabel('y', 'fontsize', 14)
zlabel('z', 'fontsize', 14)
axis equal
grid on
rotate3d on

% Add meshes legend
function lg = put_legend(hdl, names) 

    % Plot legend 
    lg = legend(hdl, names,'location','eastoutside');   
    
    % Set new position
    pos = [0.773 0.72 0.20 0.098]; 
    set(lg, 'position', pos,...
        'interpreter', 'none', 'fontsize', 11,...
        'EdgeColor', 'none', 'color', [0.98 0.96 0.96], 'autoupdate', 'off')  
    
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

% Default mesh names 
function names = default_names(N)
names = cell(N, 1);
nt = length(num2str(N));
for i = 1 : N
    sn = num2str(i);
    snum = [repmat('0', 1, nt-length(sn)), sn];
    names{i} = ['mesh-', snum];
end

function meshp = prep_patch(mesh)
% renamefields = fieldtrip function
if isfield(mesh, 'bnd')
    mesh = mesh.bnd;
end
if ~isfield(mesh, 'faces') && isfield(mesh, 'tri')
    mesh = renamefields(mesh, 'tri', 'faces');
end

if ~isfield(mesh, 'vertices')
    if isfield(mesh, 'pnt')
        mesh = renamefields(mesh, 'pnt', 'vertices');
    end
    if isfield(mesh, 'pos')
        mesh = renamefields(mesh, 'pos', 'vertices');
    end
end
meshp = [];
meshp.faces = mesh.faces;
meshp.vertices = mesh.vertices;
