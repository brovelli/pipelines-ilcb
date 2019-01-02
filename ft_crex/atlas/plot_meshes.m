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
%   - color or edgecolor : color to display the mesh edges [ 1 x 3 vector ]
%   - facecolor: specific face color [ default: [0.95 0.90 0.80] ]
%   - facealpha
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
dopt = struct('newfig', 1,...
            'dispnorm', 0,...
            'colors', [],...
            'facecolor', [],...
            'facealpha', [],...
            'edgecolor', [],...
            'edgealpha', [],... 
            'names', [],...
            'visible', 'on');

if nargin < 2
    opt = dopt;
else
    opt = check_opt(opt, dopt);
end

if ~isempty(opt.colors) && isempty(opt.edgecolor)
    opt.edgecolor = opt.colors;
end

if ischar(opt.names)
    opt.names = {opt.names};
end
opt.name = opt.names;

% Case for 'none'
if ischar(opt.facecolor)
    opt.facecolor = {opt.facecolor};
end
if ischar(opt.edgecolor)
    opt.edgecolor = {opt.edgecolor};
end
% Default appearence for meshes: colors, aplphas, names
mopt = [];
mopt.facealpha = repmat(0.15, Nm, 1);
mopt.facecolor = repmat([0.95 0.90 0.80], Nm, 1);
mopt.edgecolor = color_group(Nm);
mopt.edgealpha = repmat(0.2, Nm, 1);
mopt.name = default_names(Nm);
meshes = set_colors(meshes, opt, mopt);

if opt.newfig
    figure
    set(gcf, 'visible', opt.visible, 'units', 'centimeters', 'Position', [2 2 21 23])
    set(gca, 'position', [0.079 0.104 0.674 0.7335])
end

hold on
hnrm = struct;
hpp = zeros(Nm, 1);
names = cell(Nm, 1);
g = 1;
for i = 1 : Nm
    % Prepare mesh: identify vertices and faces fields, set color and name 
    mesh = prep_patch(meshes{i});
    if ~isempty(mesh.faces) && ~isempty(mesh.vertices)
        hpp(g) = patch('faces', mesh.faces, 'vertices', mesh.vertices,...
            'edgecolor', mesh.edgecolor, 'edgealpha', mesh.edgealpha, 'facecolor', mesh.facecolor,...
            'facealpha', mesh.facealpha, 'facelighting','gouraud', 'pickableparts', 'none'); 
        if opt.dispnorm 
            hn = add_norm(mesh.vertices, mesh.norm, mesh.edgecolor);
            if i>1
                hnrm(i) = hn;
            else
                hnrm = hn;
            end
                    
        end
        names{g} = mesh.name;
        g = g + 1;
    end   
end

hlg = put_legend(hpp(1:g-1), names(1:g-1));

xlabel('x', 'fontsize', 14)
ylabel('y', 'fontsize', 14)
zlabel('z', 'fontsize', 14)
axis equal
grid on
rotate3d on
if nargout>=1
    hp = hpp;
    hleg = hlg;
    hnorm = hnrm;
end


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
function mesh = prep_patch(mesh)

% renamefields = fieldtrip function
if isfield(mesh, 'bnd')
    bmesh = mesh.bnd;
    if isstruct(bmesh)
        fd = fieldnames(bmesh);
        Nf = length(fd);
        for i = 1 : Nf
            mesh.(fd{i}) = bmesh.(fd{i});
        end
    end
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

% Possible field name for normals
fori = {'normals', 'mom', 'ori', 'norm'};
fnam = fieldnames(mesh);
isori = ismember(fori, fnam);
if any(isori)
    mesh.norm = mesh.(fori{find(isori==1, 1, 'first')});
end
    
function meshes = set_colors(meshes, opt, mopt)
% Check each values - priority to opt parameters (even if mesh is holding the
% color field)
Nm = length(meshes);

[onames, Nc] = get_names(mopt);

for i = 1 : Nm
    msh = meshes{i};
    for j = 1 : Nc
        onam = onames{j};
        if ~isempty(opt.(onam))
            ov = opt.(onam);
            
            if length(ov)==Nm || length(ov(:, 1))==Nm
                if iscell(ov)
                    msh.(onam) = ov{i};
                else
                    msh.(onam) = ov(i, :);
                end
            else
                if iscell(ov)
                    msh.(onam) = ov{1};
                else
                    msh.(onam) = ov;
                end
            end
        else
            if ~isfield(msh, onam)
                oval = mopt.(onam);
                if iscell(oval)
                    ovals = oval{i};
                else
                    ovals = oval(i, :);
                end
                msh.(onam) = ovals;
            end
        end
    end
    meshes{i} = msh;
end
        
    

