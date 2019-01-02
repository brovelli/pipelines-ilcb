function ftData = cmeg_artefact_rm(ftData, winart)
% Removing of energetic artefacts inside continuous data
% - artefact that are visible on the major part of the channels
%
% For each channel, the artefact portions are substituted by a
% miror image of the preceding data.
% -- ftData : Fieldtrip data structure of continuous data set
% -- winart : temporal windows that define all data portions containing artefacts 
%         -> matrix [Nart x 2] col.1 : artefact window begining ; col.2 : ending
%            as many lines as artefacts to delete (Nart)
%
% If data portion is to long to be replaced by the miror image of 
% subsequent data, artefact portion values are replaced by zeros.
%
% --- CREx 2017
% CREx-BLRI-AMU project

if isempty(winart)
    return;
end

% Time vector
td = ftData.time{1};

% Number of time samples
Ns = length(td);

% Amplitudes (Nchan x time)
xall = ftData.trial{1};

% Number of artefacts
sz = size(winart);
Na = sz(1);
iart = zeros(sz);

for i = 1 : Na
    wi = winart(i, 1);
    wf = winart(i, 2);
    fprintf('\nPadding data part : [ %3.2f - %3.2f ] s\n', wi, wf);
    
    % Sample indices that fall inside artefact windows
    ipad = find(td >= wi & td <= wf);
    
    ibeg = ipad(1);
    iend = ipad(end);
    
    Nw = length(ipad);
    
    % If enough data before the artefact:
    % fill the artefact window with the mirror image of the preceding data
    if ibeg-Nw > 0
        xall(:, ipad) = xall(:, ibeg-1 : -1 : ibeg-Nw);
    else
        % Fill with the mirror image of the following data
        if iend + Nw < Ns
            xall(:, ipad) = xall(:, iend+1 : iend+Nw);
        else
            % Fill with zeros as not enough consecutive data 
            xall(:, ipad) = 0;
        end
    end
    iart(i,:) = [ibeg iend];
end

ftData.trial{1} = xall;
ftData.hdr.artpad = iart;
