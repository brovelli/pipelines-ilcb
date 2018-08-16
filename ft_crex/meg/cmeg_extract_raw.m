function rawData = cmeg_extract_raw(dpath, trl)
% Find and extract raw data found in a specific directory. 
% Raw data can be MEG 4D data or S/EEG data.
% Input paramaters:
% --------
% - dpath : path of the directory where data are looked for
%               or path of the raw data file to extract
% - trl : vector to specified trial to extract [ default: []]
%   If trl is unset, empty or 0, continuous data are extract [ default ].
%   If trl is defined (as a Nx3 matrix), trials are defined (see
%   ft_preprocessing for trl definition).
% Output parameters:
% ---------
% rawData : FieldTrip data structure containing extracted data
%
% This function uses the function filepath4D of the ft_CREx toolbox to find the 
% 4D file in the directory. Data are then extracted by ft_definetrial and
% ft_preprocessing FieldTrip functions.
%
%______
%-CREx-20180426 


%-- Check for inputs
if nargin == 1 || isempty(trl) || trl==0
    trlfield = 0;
else
    trlfield = 1;
end

% -- Search file
fprintf('\n\t-------\nSearch for raw data file\n--> %s\t-------\n', dpath);
% If dirpath is not the path of the raw file but of a directory, raw data
% are searched inside the dirpath folder using filepath4d function.
if exist(dpath,'file')==2
    datpath = dpath;
    dpath = fileparts(datpath);
else
    datpath = filepath_raw(dpath);
end

%-- Extract data with fieldtrip functions
if isempty(datpath)
    rawData = [];
    warning('Data not found in directory %s', dpath)
    return;
end

fprintf('\n\t-------\nExtract raw dataset\n\t-------\n')
[~, ~, ext] = fileparts(datpath);
% Defined type of channel to extract
% If data extent is 'eeg', then EEG channels will be extracted.
% Otherwise, data are assumed to be MEG data.
if strcmp(ext, '.eeg')
    chan = {'*'};
    styp = 'seeg';
elseif strcmp(ext, '.fif')
    chan = {'meg'}; %{'meg*1'};
    styp = 'meg';
else
    chan = {'meg'};
    styp = 'meg';
end

cfg = [];
cfg.dataset = datpath;
cfg.trialfun = 'ft_trialfun_general'; % Default (avoid warning message)
cfg.trialdef.triallength = Inf;       % All the dataset blocksize
cfg_rawData = ft_definetrial(cfg);    
cfg_rawData.channel = chan;

% If trl is define as input parameter, then data will be epoched
% according to the trl matrix.
if trlfield==1
    cfg_rawData.trl = trl;
end

rawData = ft_preprocessing(cfg_rawData);
rawData.hdr.info.rec = styp;
rawData.hdr.info.type = 'continuous_raw';

% Save the rawData mat in dpath directory if rawData is not output
if ~nargout
    save([dpath, filesep, 'rawData.mat'], 'rawData', '-v7.3')
end
