function [ftData, zchan] = cmeg_cleanup_fig(dpath, opt)
% Figures to help for bad channel identification
%
% A high pass filter with fc = 0.5 Hz is applied to the dataset to facilitate the visualization
% and the identification of recording issues that could be hidden by the low-frequency component.
%
% Figures are saved in _preproc_fig directory. 3 types of figures are
% generated:
%
% - fft plots - spectra are computed for all channels and superimposed on the
% same figure. 
%   [saved in 'fftplots' subdirectory]
%
% - continuous time series plots + 2 details of 10 s duration (one around the
% maximum amplitude). The data set is first resampled at 400 Hz to save time when recording figures
%   ['datadisp_filt' subfolder]
%
% - the layout of MEG sensor with a color code at each channel that is an
% indicator of the signal amplitude compared to the channel that displays the
% minimum value (except for the channels with zeros values)
%   ['datadisp_filt' subfolder - with the same color code being used for the time
%   series plots]
%
%
%-CREx-180726

dopt = struct('savepath', [], 'figpath', [], 'info', []);
dopt.datatyp = {'filt', 'raw'};
if nargin < 2
    opt = dopt;
else
    opt = check_opt(opt, dopt);
end

% check if dpath is a string or a data structure
if isstruct(dpath)
    ftData = dpath;
else
    [ftData, pdat] = load_data(dpath, opt);
    if isempty(opt.savepath)
        opt.savepath = dpath;
    end
    if isempty(opt.info)
        opt.info = pdat;
    end
end
    
if isempty(ftData)
    warning('Data not found to display cleanup fig\n%s\n', dpath);
end

if isfield(ftData.hdr, 'info') && (~isempty(strfind(ftData.hdr.info.type, 'raw')) ||...
        (~isempty(strfind(ftData.hdr.info.type, 'filt')) &&...
        isempty(strfind(strjoint(ftData.hdr.info.preproc.sfilt, ''), 'hp')))) %#ok
    opt = [];
    opt.type = 'hp';
    opt.fc = 0.5; 
    ftData = cmeg_filt(ftData, opt);
end

spath = opt.savepath;

if isempty(opt.figpath)
    pfig = make_dir([spath, filesep, '_preproc']);
else
    pfig = opt.figpath;
end

% Subdirectory for fft plots
pfig_sp = make_dir([pfig, filesep,'fftplots'], 1);    

% Launch fft calculations
% Stacked spectrum with default paramaters (spsparam.n=30 and
% spsparam.dur=20 s) :
spData = meg_fftstack_calc(ftData);    

% Figures of superimposed stacked spectrum per channels
if ~isempty(spData.spectra)
    param = spData.param;
    % Figures with the dedicated subfunction
    fopt = opt;
    fopt.savepath = pfig_sp;
    fopt.sps_param = param;
    meg_fftstack_fig(spData, fopt) 
end

% Figure avec des zooms dans les donnees pour mieux voir si probleme
pfig_vis = make_dir([pfig, filesep, 'datadisp'], 1); 

Fs = fsample(ftData.time{1});
if Fs > 400
    cfg = [];
    cfg.resamplefs = 400;
    ftData = ft_resampledata(cfg, ftData);    
end

[Scol, zchan] = cmeg_chancheck_fig(ftData, pfig_vis, opt.info);

fopt = [];
fopt.savepath = pfig_vis;
fopt.info = opt.info;
fopt.colors = Scol;
cmeg_plot_data(ftData, fopt)


function [ftData, pdat] = load_data(dpath, opt)
% Find data accoprding to opt.datatyp
pdat = find_mat(dpath, opt);
if isempty(pdat)
    pdat = [dpath, filesep, 'rawData'];
    % Try to extract from raw data
    ftData = cmeg_exctract_raw(dpath);
else    
    ftData = loadvar(pdat);
end

