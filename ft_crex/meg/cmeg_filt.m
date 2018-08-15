function [filtData, fopt] = cmeg_filt(rawData, fopt)
% Sub-function to launch the filtering by ft_preprocessing FieldTrip function
% The order of the filters are included as fixed parameters with a value = 2.
% If no cut-off frequency is given in the filtopt.fc field or if only one
% value is store for band-pass filter, the missing values are asking on the 
% command window.
% 
%-CREx 20140520 
%-CREx-BLRI-AMU project: https://github.com/blri/CREx_MEG_FieldTrip


if nargin < 1 
    fopt = [];   
end

if ~isempty(fopt) && isfield(fopt, 'channel')
    chan = fopt.channel;
else
    chan = '*';
end

fopt = cmeg_filt_opt(fopt);

fc = fopt.fc;
fopt = fopt.type;

cfg = [];
cfg.channel = chan;
cfg.instabilityfix = 'split';
switch lower(fopt)
    case 'none'
        if ~iscell(chan)
            filtData = rawData;
            return
        end
        
    case 'hp'
        cfg.hpfilter = 'yes';
        cfg.hpfiltord = 2;
        cfg.hpfreq   = fc;
        
    case 'lp' 
        cfg.lpfilter = 'yes';
        cfg.lpfiltord = 2;
        cfg.lpfreq   = fc;
        
    case 'bp'  
        cfg.bpfilter = 'yes';
        cfg.bpfiltord = 2;
        cfg.bpfreq   = fc;    
end

filtData = ft_preprocessing(cfg, rawData);