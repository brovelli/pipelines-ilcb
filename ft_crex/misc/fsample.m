function Fs = fsample(vtime)
% Compute the sampling frequency from a time vector vtime
% If vtime is a structure, the value of the field 'time' is set as time vector.
% If vtime is a cell, only the first element (vtime{1}) is used.
%
%-CREx-180425

if isstruct(vtime)
    vtime = get_field(vtime, 'time');
end

if iscell(vtime)
    vtime = vtime{1};
end
    
Nt = length(vtime);
if Nt > 1
    Fs = (Nt - 1) / diff(vtime([1 end]));
else
    Fs = [];
    warning('Sampling frequency not determinable');
end

