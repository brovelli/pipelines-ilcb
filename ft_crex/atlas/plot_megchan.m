function plot_megchan(grad)

% Be sure units in mm to match with fwd files
grad = ft_convert_units(grad, 'mm');
    
% Select only MEG label A* for plotting
clab = char(grad.label);
idplot = strfind(clab(:,1)','A');

% Set coordinate units to mm (from m)
chanpos = grad.chanpos; 
chanpos = chanpos(idplot,:);
p = plot3(chanpos(:,1),chanpos(:,2), chanpos(:,3),'o');
set(p,'markersize', 6,'markerfacecolor', [0 0.6 0.75], 'markeredgecolor', [1 0.7 0])