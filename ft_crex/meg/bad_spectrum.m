function zchan = bad_spectrum(spData, ftData, opt)
dopt = [];
dopt.thr = 4;
% Frequency band to check for mean spectra value
dopt.bd = [2 200];

if nargin < 3
    opt = dopt;
else
    opt = check_opt(opt, dopt);
end

% Identify spectrum with extremum value in fbd frequency band
m = opt.thr;
fbd = opt.bd;

freq = spData.freq;
spectra = spData.spectra;

% Mean spectrum values in frequency band
msp = mean(spectra(:, freq >fbd(1) & freq < fbd(2)), 2);
% Standard deviation of mean values
sdev = std(msp);

% Index of extremum values (median value +/- a factor (m) * sdev)
ibad = find(msp <= median(msp)-m*sdev | msp > median(msp)+m*sdev);
if isempty(ibad)
    zchan = [];
else
    zchan = spData.label(ibad);
end

% Identify channels with too much 0
xd = ftData.trial{1};
isz = sum(xd==0, 2) > length(xd(1,:))/2;
if any(isz)
    zchan = unique([zchan; ftData.label(isz)]);
end