function pdat = filepath_raw(dirpath, dtyp)
% Return full access path to MEG (4D data) or S/EEG data in the
% directory specified by 'dirpath' path.
% dtyp is an optionnal parameter that indicates the data type to find in the
% search directory.
%
% dtyp values according to the requested data file type:
%   'raw'    : raw data ("c,rfDC")
%   'rawcor' : raw data with noise correction ("c,rfDC,n") [default]
%   'filt'   : filtered data ("c,rf*Hz")
%   'filtcor': filetered data with noise correction ("c,rf*Hz*n")
%   'seeg'   : S/EEG data with .eeg file extent
%
% If dtyp string is not supplied or empty, data will be searched with this
% default order (as soon as corresponding data is found, searching stops):
% (1) 'c,rfDC,n' : MEG 4D raw data with noise reduction (<=> dtyp = 'rawcor') 
% (2) 'c,rfDC'   : MEG 4D raw data ('raw')
% (3) 'c,rf*Hz*n': MEG 4D filtered with noise reduction ('filtcor')
% (4) 'c,rf*Hz'  : MEG 4D filtered ('filt')
% (5) '*.eeg'    : S/EEG data ('seeg')
% (6) 'fif'      : MEG Neuroimagery
% The path of the first found file is return.
%
%______
%-CREx 20131030 

fprintf('\n-----\nLooking for raw data set in %s\n', dirpath);

% List of possible data set
raw_nam = {'c,rfDC,n', 'c,rfDC', 'c,rf*Hz*n', 'c,rf*Hz', '*.eeg', '*.fif'}; 
% Associated data type
dtyp_nam = {'rawcor', 'raw', 'filtcor', 'filt', 'seeg', 'rawns'};

pdat = [];
% If dtyp is specified, looking for the corresponing raw data set
if nargin==2 && ~isempty(dtyp) && any(strcmp(dtyp, dtyp_nam))
    datnam = raw_nam{strcmp(dtyp, dtyp_nam)==1};
    dpath = fullfile(dirpath, datnam);
    dp = dir(dpath);
    if isempty(dp)
        warning('!!');
        fprintf('4D (or seeg) data "%s" not found in directory\n%s\n',datnam, dirpath);
    else
        pdat = fullfile(dirpath, dp(1).name);
    end
else
    % Default raw data set are reseached in dirpath directory
    Nr = length(raw_nam);
    for i = 1 : Nr
        dpath = fullfile(dirpath, raw_nam{i});
        dp = dir(dpath);        
        if ~isempty(dp)
            pdat = fullfile(dirpath, dp(1).name);
            break;
        end
    end
end
if ~isempty(pdat)
    [~,fnam] = fileparts(pdat); 
    fprintf('\n---- Found : %s\n\n',fnam);
else
    warning('Raw data set not find in directory')
    fprintf('\nCheck for the raw data list inside filepath_raw script (raw_list, line 30)\n!!!!\n\n');
end