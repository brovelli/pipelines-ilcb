function [layout, unit] = cmeg_det_lay(label)
%- Try to determine the fieldtrip layout to use for plotting from the
% sensor label name

slab = label{1};
if length(slab)> 3 && strcmp(slab(1:3), 'MEG')
    layout = 'neuromag306mag.lay';
    unit = 'T';
else
    ucc = unique(cellfun(@(x) x(1), label));
    Nu = length(ucc);
    % 4D MEG data with label beginning by 'A'
    if Nu==1 && length(label) > 220
        layout = '4D248_helmet.mat';
        unit = 'T';
    else
        % Probably SEEG data
        layout = [];
        unit = 'V';
    end
end