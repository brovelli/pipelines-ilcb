function Sbals = dipole_spheres(dip_pos, lab, radius, col, nfaces, alpha)
% Define the sphere patches at dipole locations 
% dip_pos : 3D dipole coordinates as Ndip x 3 matrix
%
% To plot the dipoles spheres: 
% figure, hold on
% arrayfun(@patch, Sbals)
%
%- CREx-181002

if ~nargin || isempty(dip_pos)
    Sbals = [];
    warning('No dipole coordinates provided to define the sphere patchs')
    return
end
% Number of faces to define the dipole sphere
if nargin < 5 || isempty(nfaces)
    nfaces = 10;
end
if nargin < 6 || isempty(alpha)
    alpha = 1;
end
% Number of dipoles 
Ndip = size(dip_pos, 1);

% Color for each dipole [default: all in red]
% A Ndip x 3 RGB color (with values from 0 to 1) or a 1 x 3 to apply same color
% for all dipoles
if nargin < 4 || isempty(col)
    col = [0.85 0 0];
end
    
if Ndip > 1 && length(col(:, 1))==1 
    col = repmat(col, Ndip, 1);
elseif length(col(:, 1)) ~= Ndip
    warning('Color matrix size not consistent with dip_pos size (Ndip x 3)')
    warning('All dipole patch will be defined with the same color')
    col = repmat(col(1, :), Ndip, 1);
end

if nargin < 3 || isempty(radius)
    radius = 1;
end

if Ndip > 1 && length(radius)==1 
    radius = repmat(radius, Ndip, 1);
elseif length(radius) ~= Ndip
    warning('Color matrix size not consistent with dip_pos size (Ndip x 3)')
    warning('All dipole patch will be defined with the same color')
    radius = repmat(radius(1), Ndip, 1);
end

if nargin < 2 || isempty(lab) || length(lab)~=Ndip
    lab = repmat({'dip_'}, Ndip, 1);  
    lab = cellfun(@(x,y) [x, num2str(y)], lab, num2cell(1:Ndip)', 'UniformOutput', false);
end

%----
% Set dipoles spheres
[Xsph, Ysph, Zsph] = sphere(nfaces); 
% Unit sphere with N=nfaces faces - Xb, Yb and Zb are (N+1) x (N+1) matrices
% using by surf(Xb, Yb, Zb) to draw the sphere
% Set radius of sphere representing dipoles = 6 mm

%----
% Define dipoles patch balls properties

Sbals = struct;
for i = 1: Ndip 

    % Coordinates of the center
    cx = Xsph.*radius(i) + dip_pos(i, 1);
    cy = Ysph.*radius(i) + dip_pos(i, 2);
    cz = Zsph.*radius(i) + dip_pos(i, 3);

    [Sbals(i).faces, Sbals(i).vertices] = surf2patch(cx, cy , cz );   
    Sbals(i).facecolor = col(i, :);
    Sbals(i).edgecolor = 'none';
    Sbals(i).facelighting = 'gouraud';
    Sbals(i).displayname = lab{i};
    Sbals(i).facealpha = alpha;
end