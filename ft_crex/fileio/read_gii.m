function dat = read_gii(pdat)
% Read gii data for pipline project
%-CREx-20180206

% Add gifti toolbox (fieldtrip/external/gifti)
ft_tool('gifti', 'add');

gdata = gifti(pdat);

fnames = fieldnames(gdata);

% Case mesh
if any(strcmp(fnames, 'vertices'))
    dat = export(gdata, 'ft');
end

% Case texture
if any(strcmp(fnames, 'cdata'))
    dat = double(gdata.cdata);
end

% Remove gifti toolbox 
% ft_tool('gifti', 'rm');