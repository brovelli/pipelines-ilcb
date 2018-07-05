function labelinfo = read_labelinfo(pinfo)
% Read the labelinfo.txt file containing MarsAtlas label information
% Parse the Brodman area
%
% pinfo: path to the labelinfo.txt file [default: those in ft_crex/atlas toolbox]
%
% labelinfo: structure with fields
%   - index: index of the atlas parcel
%   - hemis: hemisphere (L or R)
%   - lobe: lobe name (ex. Frontal)
%   - code: code name (the same for the 2 hemisphere) (ex. Mdl)
%   - fullname (ex. DorsolateralMotorCortex)
%   - label: code name with hemisphere indication (ex. Mdl_L)
%   - BA: Brodman area n° (ex. 4) (empty or vector of numbers)
%
% If nargout==0, then the file ft_crex/atlas/labelinfo.mat is update with the
% labelinfo structure.
% 
%
%-CREx-180702

if ~nargin || isempty(pinfo)
    pinfo = fullfile(ptool, 'atlas', 'labelinfo.txt');
end

fid = fopen(pinfo, 'r');

shdr = fgetl(fid);
Cdat = textscan(fid, '%d%s%s%s%s%s%s', 'Delimiter', '\t');
fclose(fid);

cnam = strsplitt(shdr, '\t');

idx = Cdat{strcmp(cnam, 'Index')};

labelinfo = [];
labelinfo.index = idx;
labelinfo.hemis = Cdat{strcmp(cnam, 'Hemisphere')};
labelinfo.lobe = Cdat{strcmp(cnam, 'Lobe')};
labelinfo.code = Cdat{strcmp(cnam, 'Code')};
labelinfo.fullname = Cdat{strcmp(cnam, 'FullName')};
labelinfo.label = Cdat{strcmp(cnam, 'Label')};
labelinfo.BA = Cdat{strcmp(cnam, 'BA')};

% Remove the WhiteMatter (idx==0)
labelinfo = structfun(@(x) x(idx~=0), labelinfo, 'UniformOutput', 0);

% Parse the Brodman Areas
labelinfo.BA = cellfun(@(x) parse_ba(x), labelinfo.BA, 'UniformOutput', 0);

if ~nargout
    save(fullfile(ptool, 'atlas', 'labelinfo.mat'), 'labelinfo')
end


function vba = parse_ba(sba)
vba = [];
if ~strcmp(sba, 'NA')
    cba = strsplitt(sba, '/');
    vba = cellfun(@(x) str2double(x), cba);
end
        