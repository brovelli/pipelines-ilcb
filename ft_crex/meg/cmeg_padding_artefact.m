function padData = cmeg_padding_artefact(ftData, windef)
% Removing of energetic artefacts inside continuous data
% - artefact that are visible on the major part of the channels
% For data at each channels, the artefact portions are substituted by a
% miror image of the preceding data. All the events falling inside the 
% artefact portions are removal from the cfg_event structure (that will be
% use to do data epoching by meg_ft_preproc_extractrial)
% -- ftData : Fieldtrip data structure of continuous data set
% -- windef : temporal windows that define all data portions containing artefacts 
%         -> matrix [Nart x 2] col.1 : artefact window begining ; col.2 : ending
%            as many lines as artefacts to delete (Nart)
%
% If data protion is to long to be replaced by the miror image of 
% subsequent data, artefact portion values are replace by zeros.
% --- CREx 2017
% CREx-BLRI-AMU project

if isempty(windef)
    padData = ftData;
    return;
end

% Time vector
td = ftData.time{1};

% Number of time samples
Ns = length(td);

% Amplitudes (Nchan x time)
xall = ftData.trial{1};

% Number of artefacts
Na = length(windef(:,1));
artpad = zeros(size(windef));

for i = 1 : Na
    wi = windef(i, 1);
    wf = windef(i, 2);
    fprintf('\nPadding data part : [ %f - %f ] s\n', wi, wf);
    % Sample indices that fall inside artefact windows
    ipad = find(td >= wi & td <= wf);
    
    ibeg = ipad(1);
    iend = ipad(end);
    
    Nwp = length(ipad);
    
    % If enough data before the artefact:
    % fill the artefact window with the mirror image of the preceding data
    if ibeg-Nwp > 0
        xall(:,ipad) = xall(:, ibeg-1 : -1 : ibeg-Nwp);
    else
        % Fill with the mirror image of the following data
        if iend + Nwp < Ns
            xall(:,ipad) = xall(:, iend+1 : iend+Nwp);
        else
            % Fill with zeros as not enough consecutive data 
            fprintf('\n--- !!\nNot enough duration before the artefact to pad\n');
            fprintf('artefact with previous mirror data image : \n');
            fprintf('Padding artefact window with ZEROS\n');
            xall(:,ipad) = 0;
        end
    end
    artpad(i,:) = [ibeg iend];
end

padData = ftData;
padData.trial{1} = xall;

padData = add_info(padData, 'pad', artpad);
padData.hdr.info.artpad = artpad;