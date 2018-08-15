function spData = meg_fftstack_calc(ftData, sps_param)
% Calcul d'un spectre moyen pour chaque capteur 
%
% Un spectre moyen ("stacked spectrum") corresponda la somme de plusieurs 
% spectres calcules sur des portions de donnees consecutives de longueur fixe, 
% divisee par le nombre de portion. 
%
%____ INPUT
%
% ftData : data structure from FieldTrip processing of continuous data
% sps_param : parameters for the stacked spectra calculation
%   sps_param.dur : data duration for each elementary  pour le calcul de chaque spectre
%      elementaire, en secondes
%   spsparam.n   : nombre de spectres elementaires a calculer sur les portions 
%      consecutives du signal de duree spsparam.dur
%
%____ OUTPUT
%
% allFFTstack : les spectres moyens pour tous les capteurs (un spectre par
%   ligne = pour un capteur)
% freqst : le vecteur des frequences associe aux valeurs d'amplitude des
%   spectres contenus dans allFFTstack
%
% Trace des spectres moyens de tous les capteurs :
% >> plot(freqst,allFFTstack)
% Trace du spectre moyen calcule sur les donnees du capteur
% FTData.label{10} :
% >> plot(freqst,allFFTstack(10,:))
% 
% Un spectre moyen represente le contenu frequentiel moyen d'une partie
% des donnees de duree spsparam.n x spsparam.dur secondes. 
%
% Valeurs par defaut si spsparam n'est pas entre en argument de la fonction :
% spsparam.dur = 20 et spsparam.n = 30
% => Un spectre stacke est obtenu par le moyennage de 30 spectres calcules 
% sur chaque portion consecutive des donnees temporelles de 20 secondes. 
% 
% Pour exclure les possibles artefacts au niveau de la bordure inferieur des
% donnees, le calcul des fft ne debute qu'apres les 2 premieres portions de
% spsparam.dur secondes.
% => avec les parametres par defaut par exemple, les calculs de spectres 
% est effectue entre 40 s et 340 s.
%
% Si les parametres spsparam entres par l'utilisateur ne conviennent pas,
% les parametres par defaut sont testes. S'il ne conviennent pas a leur
% tour, le programme retourne des valeurs vides de feqst et allFFTstack.
% (ex. : duree des donnees trop courte pour calculer spsparam.n spectres
% sur des portions de spsparam.dur s)
% 
%
% Les spectres sont calcules par l'algorithme FFT (Fast Fourier Transform), 
% implemente dans la fonction fft de Matlab (Signal Processing Toolbox). 
%
%----------
% CREx 20131129

% Initialized
spData = [];
spData.spectra = [];
spData.freq = [];
if ~isfield(ftData, 'fsample')
    fsamp = fsample(ftData.time{1});
else
    fsamp = ftData.fsample;
end

defp = [];
% Elementary duration in second
% Duration of each data portion to compute the spectrum
defp.dur = 20;
% Number of elementary and consecutive spectum to average
defp.n = 30;

% Apply default parameters if sps_param is not set
if nargin < 2 || isempty(sps_param) || ~isstruct(sps_param)
    sps_param = defp;  
else
    % Check for the sps_param structure
    if ~isfield(sps_param, 'dur') || ~isnumeric(sps_param.dur) || length(sps_param.dur)~=1
        % Default parameter for duration
        sps_param.dur = defp.dur;
    end
        
    if ~isfield(sps_param,'n')  || ~isnumeric(sps_param.n) || length(sps_param.n)~=1 
        % Default parameter for number
        sps_param.n = defp.n;   
    end
end 

xall = ftData.trial{1};
Nchan = length(xall(:,1));
Ndat = length(xall(1,:));

% Check for parameters according to data length
sps_param = check_stackparam(Ndat, fsamp, sps_param, defp);

if sps_param.n < 2
    fprintf('\n!!!\nNumber of consecutive spectrum to average is too small (< 2)\n')
    fprintf('Calculation is cancelled\n')
    return
end

idec = sps_param.int;
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
spData.param = sps_param;

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
fprintf('    %s s-duration consecutive portions\n    of data will be stacked', num2str(param.dur));
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

    

