function masurf = read_mars_surf(psurf, ptex)
% Read surfaces and associated textures from BV / MarsAtlas processing

if nargin < 2
    ptex = [];
else
    if ~iscell(ptex)
        ptex = {ptex};
    end
end

if ~iscell(psurf)
    psurf = {psurf};
end

% Number of surface (+ associated textures) to read
Nh = length(psurf);

masurf = cell(Nh, 1);
for k = 1 : Nh
    psurfh = psurf{k};
    [~, snam] = fileparts(psurfh);
    % Define if Left or Right according to surf name ('Lwhite' or 'Rwhite')
    % Normally, Lwhite is in k==1 and Rwhite in k==2 as path list is
    % returned in alphabetic order - but in case of...
    if ~isempty(strfind(snam, 'Lwhite')) %#ok
        % snam = 'Lwhite';
        namh = 'surf_L';
    else
        % snam = 'Rwhite';
        namh = 'surf_R';
    end
    masurfh = read_gii(psurfh);
    % ptex list of file in the same order as psurf a priori (if not, see to
    % find the corresponding ptex depending on snam)
    if ~isempty(ptex)
        masurfh.tex = read_gii(ptex{k}); 
    end
    masurfh.name = namh;
    masurf{k} = masurfh;        
end