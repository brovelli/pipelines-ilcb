function dat = read_gii(pdat)
% Read gii data for pipeline project (GII surface and texture files)
%
% Use gifti function included in fieldtrip/external/gifti toolbox to extract data from gifti object
% 
% For mesh data (surf*.gii files), the normals are added to the dat structure if
% it was not extracted from the raw GII file, using pathnormals function 
% (see ft_crex/external/patch_normals).
%
%-CREx-180206

% Add gifti toolbox (fieldtrip/external/gifti)
ft_tool('gifti', 'add');

gdata = gifti(pdat);

% Case mesh
if isfield(gdata, 'vertices')
    dat = [];
    dat.pos = double(gdata.vertices);
    dat.tri = double(gdata.faces);

    % Extract mat field leads to an error
    %     Error in gifti/subsref (line 28)
    %                 varargout{1} = this.data{j}.space.MatrixData;
    %     if isfield(gdata, 'mat')
    %         dat.mat = double(gdata.mat);
    %     end
    
    % To check if it doesn't bug when normals are actually included in gii file
    %     if isfield(gdata, 'normals')
    %         dat.ori = double(gdata.normals);
    %     end
else
    % Case texture
    if isfield(gdata, 'cdata')
        dat = double(gdata.cdata);
    end
end

% Remove gifti toolbox 
% Important to unable the 'isfield' function of gifti toolbox
ft_tool('gifti', 'rm');