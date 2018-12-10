function ofreq = prep_freq_param(freq_param, pres)
% Frequency analysis parameters for ft_freqanalysis from cp_main bopt parameters
%-CREx180920

cfg = [];
cfg.method = 'mtmconvol';
cfg.output = 'fourier';
cfg.taper = 'dpss';
cfg.keeptrials = 'yes';

Nf = length(freq_param(:, 1));

for i = 1 : Nf
    fparam = freq_param(i, :);
    foi = fparam{2};
    nbd = length(foi);
    Df = fparam{3};
    Dt = fparam{4};
    if nbd > 1 && length(Df)==1
        Df = ones(1, nbd)*Df;
    end
    if nbd > 1 && length(Dt)==1
        Dt = ones(1, nbd)*Dt;
    end
    
    % Frequency analysis directory
    [dana, dfreq] = dir_freq(fparam);
    
    pfrq = make_dir([pres, filesep, dana]);
    % More specific directory
    pfrq = make_dir([pfrq, filesep, dfreq]);
    
    % Number of Slepian tapers (3 is good)
    ntap = 2 .* Dt .* Df - 1;  
        
    So = cfg;
    
	%- frequency(ies) of interest in Hz
    So.foi = foi;
    %- width of the frequency smoothing in Hz = frequency resolution Df 
    So.tapsmofrq = Df;
    %- time window width (time resolution) in s
    So.t_ftimwin = Dt;
    %- times of interest in s = middle time locations of the sliding window 
    So.toi = [];
           
    % Additionnal info
    
    % Directory for results
    So.dir = pfrq;
    
    % Keep param for result mat
    param = [];
    % Name of freq analysis
    param.name = dana;
    param.df = Df(1);
    param.dt = Dt(1);
    param.foi = foi;
    param.ntap = ntap;
    
    So.param = param;
    if i==1
        ofreq = So;
    else
        ofreq(i) = So;
    end
end

% Directory names according to frequency parameters
function [dana, dfreq] = dir_freq(freqpar)
% Freq analysis directory name (ex. 'hga')
dana = name_save(freqpar{1});
% Specific freq parameter directory name
foi = num2str(freqpar{2});
df = num2str(freqpar{3});
dt = strrep(num2str(freqpar{4}, '%1.3f'), '.', '');
dfreq = strrep(['fq', foi, '_df', df, '_dt', dt], '.', 'p');
