function spData = cmeg_fft_mean(ftData, opt)
% Mean spectrum calculation on continuous data for each sensor
%
% A mean spectrum = mean of spectra computed from consecutive portions of data 
%____ INPUT
%
% ftData : data structure from FieldTrip processing of continuous data
% opt : parameters for the mean spectra calculation
%   opt.dur : data duration for each elementary spectrum calculation in second
%   opt.n   : number of consecutive elementary spectrum for averaging
%
%----------
% CREx 20131129

% Initialized
spData = [];

if ~isfield(ftData, 'fsample')
    fsamp = fsample(ftData.time{1});
else
    fsamp = ftData.fsample;
end

dopt = [];
% Elementary duration in second
% Duration of each data portion to compute the spectrum
dopt.dur = 20;
% Number of elementary and consecutive spectum to average
dopt.n = 30;

if nargin < 2
    opt = dopt;
else
    opt = check_opt(opt, dopt);
end

xall = ftData.trial{1};
Nchan = length(xall(:,1));
Ndat = length(xall(1,:));

% Check for parameters according to data length
opt = check_stackparam(Ndat, fsamp, opt, dopt);

if opt.n < 2
    fprintf('\n!!!\nNumber of consecutive spectrum to average is too small (< 2)\n')
    fprintf('Calculation is cancelled\n')
    return
end

idec = opt.int;
% Number of point to compute fft
nsamp = 2^nextpow2(idec(2)-idec(1));	
% Frequencies vector
freq = (fsamp/nsamp)*(0:(nsamp/2)-1);
Nf = length(freq);

Ndec = length(idec) - 1;

spectra = zeros(Nchan, Nf);
for j = 1 : Ndec
    int = idec(j) : idec(j+1);
    % Spectra
    sp = fft(xall(:, int), nsamp, 2);
    % Normalized spectra 
    spnorm = 2*abs(sp(:, 1 : nsamp/2))./nsamp;
    spectra = spectra + spnorm;
end

spectra = spectra./Ndec;
spData.spectra = spectra;
spData.freq = freq;
spData.label = ftData.label;
spData.param = opt;

% Check for parameters for processing the stacked spectrum
% Reduce the number of spectrum to average according to the actual data length
function param = check_stackparam(Ndat, fsamp, param, defparam)

fprintf('\nCheck for stacked spectra calculation parameters\n');

intdat = fft_param(param, fsamp);

% Duration of data too short
if intdat(end) > Ndat
    fprintf('Data duration too short: reduced number of consecutive spectra to average\n');

    % Reduce the number of consecutive spectra to be stacked on regard to data
    % length
    % Number of sample per elementary spectrum is keeping
    param.n = floor((Ndat-intdat(1))./(fsamp*param.dur));
    fprintf('Set number of consecutive spectra to %f\n', param.n);
    % New vector of indices
    intdat = fft_param(param, fsamp);
    % But number of spectra to average should be > 2
    % Try with default parameters if not the case 
    if param.n < 2
        fprintf('\n\nNumber of consecutive elementary spectrum too small (< 2)\n')
        fprintf('Try by setting default parameters\n')
        fprintf('-- Duration of elementary spectrum = %f s\n', defparam.dur);
        fprintf('-- Total number of spectrum = %f \n', defparam.n);            
        intdat = fft_param(defparam);
        % Still too much spectra number with default duration 
        if intdat(end) > Ndat
            % Reduction of the number of consecutive spectrum to be stacked
            % again
            param.n = floor((Ndat-intdat(1))./(fsamp*defparam.dur));
            intdat = fft_param(param, fsamp); 
            fprintf('\n\nNumber of consecutive elementary spectrum too small (< 2)\n');
            fprintf('Set number to %f\n', param.n);
        end     
    end 
end

param.int = intdat;
fprintf('\n- - - - - - -\n => %s spectra calculated on each\n', num2str(param.n));
fprintf('    %s s-duration consecutive portions\n    of data will be averaged', num2str(param.dur));
fprintf('\n- - - - - - -\n');

% Compute intervals of indices to define consecutive portions of data to
% calculate spectra
function intdat = fft_param(param, fsamp)

% Number of sample per elementary spectrum
nusp = floor(param.dur*fsamp+1);  

% Total length of used data to make the stack
ntot = nusp*param.n;  

% The first sps_dur*2 s are excluded 
ixi = floor(param.dur*2*fsamp+1); 
ixf = ixi + ntot-1;
intdat = ixi : nusp : ixf;

    

