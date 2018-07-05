function [hp, hleg, hnorm] = plot_meshes(meshes, opt)
% Plot one or several meshes on the same figure
%
%*** Inputs
%
% -- meshes : cell of meshes or unique mesh data [REQUIRED]
% 
% Each mesh that is containing in meshes cell is a patch structure with fields:
%
%   - vertices (can be 'pos', 'pnt' or 'vertices')
%   - faces ('tri' or 'faces')
%   - if 'bnd' field is found, the vertices and faces are searched in mesh.bnd
% 
%   Optionnal fields:
%
%   - normals (can be 'ori', 'mom', 'norm' or 'normals')
%   - name: name of the mesh
%   - color: color to display the mesh [ 1 x 3 vector ]
% 
%
% -- opt: structure of options with fields
%
%   - newfig: flag to indicate to generate of new figure for plotting the meshes 
%           [ default: 1 : new figure is displayed ]
%
%   - dispnorm: flag to indicated if surface normals are to be displayed if
%           mesh contains normals (in 'mom', 'ori' or 'normals' fields)
%           [ default: 0 (do not display the normals) ]
%
%   - names : cell of mesh names for multi-mesh plot or string name for individual mesh 
%           The name of each mesh can be set as an additionnal field in each mesh
%           structure.  If opt.names is provided, the mesh.name field will not be
%           considered.
%           [default : mesh_1, mesh_2...] 
%
%   - colors : associated colors as a [Nmesh x 3] matrix 
%           In no opt.colors is set, the mesh.color field is checked to define
%           the color of the mesh
%           [default : color defined by color_group(Nmesh) function]
%
%   - visible : if newfig==1, the visible property for gcf is set to opt.visible
%       [ default: 'on']
%
% Use ft_CREx color_group function and the renamefields fieldtrip function
%
%-CREx-180704

if ~iscell(meshes)
    meshes = {meshes};
end

% Remove empty meshes
meshes = meshes(~cellfun(@isempty, meshes));
Nm = length(meshes);

% Set default
dopt = struct('newfig', 1, 'dispnorm', 0, 'colors', [], 'names', [], 'visible', 'on');
if nargin < 2
    opt = dopt;
else
    opt = check_opt(opt, dopt);
end

% Default colors
if isempty(opt.colors)
    opt.colors = color_group(Nm);
    fcol = 0;
else
    % Force color to be the one set in opt.col even if a mesh.color field is
    % found
    fcol = 1;
end

if isempty(opt.names)
    opt.names = default_names(Nm);
    fnam = 0;
else
    fnam = 1;
end

if ischar(opt.names)
    opt.names = {opt.names};
end

if opt.newfig
    figure
    set(gcf, 'visible', opt.visible, 'units', 'centimeters', 'Position', [2 2 21 23])
    set(gca, 'position', [0.079 0.104 0.674 0.7335])
end

hold on
hnorm = [];
hp = zeros(Nm, 1);
names = cell(Nm, 1);
g = 1;
for i = 1 : Nm
    % Prepare mesh: identify vertices and faces fields, set color and name 
    mesh = prep_patch(meshes{i}, opt.names{i}, fnam, opt.colors(i, :), fcol);
    if ~isempty(mesh.faces) && ~isempty(mesh.vertices)
        hp(g) = patch('faces', mesh.faces, 'vertices', mesh.vertices, 'edgecolor', mesh.color, 'edgealpha', 0.2, 'facecolor', [0.95 0.90 0.80],...
            'facealpha', 0.15, 'facelighting','gouraud', 'pickableparts', 'none'); 
        if opt.dispnorm 
            hnorm = add_norm(mesh.vertices, mesh.norm, mesh.color);    
        end
        names{g} = mesh.name;
        g = g + 1;
    end
    
end

hleg = put_legend(hp(1:g-1), names(1:g-1));

xlabel('x', 'fontsize', 14)
ylabel('y', 'fontsize', 14)
zlabel('z', 'fontsize', 14)
axis equal
grid on
rotate3d on

% Add normals
function hq = add_norm(pos, vec, meshcol)
if isempty(vec)
    return;
end
hq = quiver3(pos(:, 1), pos(:, 2), pos(:, 3), vec(:, 1), vec(:, 2), vec(:, 3));
set(hq, 'AutoScaleFactor', 1, 'linewidth', 0.7);

% Determin a darker color than mesh color
[mcol, i] = max(meshcol);
dcol = meshcol - 0.2;
dcol(i) = mcol + 0.1;
dcol(dcol < 0) = 0;
dcol(dcol > 1) = 1;

set(hq, 'Color', dcol);

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

% Prepare mesh data for plotting
function meshp = prep_patch(mesh, dname, fn, dcol, fc)
% Keep initial mesh
meshi = mesh;
% renamefields = fieldtrip function
if isfield(mesh, 'bnd')
    mesh = mesh.bnd;
end
if ~isfield(mesh, 'faces') 
    if isfield(mesh, 'tri')
        mesh = renamefields(mesh, 'tri', 'faces');
    else
        mesh.faces = [];
    end
end

if ~isfield(mesh, 'vertices')
    if isfield(mesh, 'pnt')
        mesh = renamefields(mesh, 'pnt', 'vertices');
    elseif isfield(mesh, 'pos')
        mesh = renamefields(mesh, 'pos', 'vertices');
    else
        mesh.vertices = [];
    end
end

% Set the prepare mesh structure
meshp = [];
meshp.faces = mesh.faces;
meshp.vertices = mesh.vertices;
meshp.norm = [];

% Possible field name for normals
fori = {'normals', 'mom', 'ori', 'norm'};
fnam = fieldnames(mesh);
isori = ismember(fori, fnam);
if any(isori)
    meshp.norm = mesh.(fori{find(isori==1, 1, 'first')});
end

if ~fn && isfield(meshi, 'name')
    meshp.name = meshi.name;
else
    meshp.name = dname;       
end

if ~fc && isfield(meshi, 'color')
    meshp.color = meshi.color;
else
    meshp.color = dcol;
end
    

