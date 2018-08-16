function cchan = chan_sel(crms)
% Define the channel cell for ft_channelselection from the bad channel cell
% (cfg.channel = cchan)
if isempty(crms)
    cchan = '*';
else    
    crms = cellfun(@(x) ['-', x], crms, 'UniformOutput', 0);
    cchan = ['*', crms'];
end
    
