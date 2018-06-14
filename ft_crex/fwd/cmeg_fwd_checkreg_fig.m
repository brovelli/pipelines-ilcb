function cmeg_fwd_checkreg_fig(vol, grad, infodata, savpath)
% Check coregistration figure
% vol : head model volume (singleshell)
% grad : MEG sensor grad structure with sensors coordinates
% infodata : additionnal information to add on the figure (in figure title) [default : '']
% savpath : directory path where to save the figure (if not specified, figure
%           will be generated but not saved) [default: []]
%
%-CREx-170228

if nargin < 3
    infodata = [];
end
if nargin < 4 || isempty(savpath)
    sav = 0;
else
    sav = 1;
end
% Select only MEG label A* for plotting
clab = char(grad.label);
idplot = strfind(clab(:,1)','A');

% Set coordinate units to mm (from m)
chanpos = grad.chanpos.*1000; 
chanpos = chanpos(idplot,:);

figure
set(gcf, 'units','centimeters','position',[5 6 20 18])
%-- [1] Head volume + channel in x-z plan
subplot(221)
checkreg_plot(vol, chanpos)
view(0,0)
title('From the right')
xlabel('x (mm)')
zlabel('z (mm)')

%-- [2] Head volume + channel in y-z plan
subplot(222)
checkreg_plot(vol, chanpos)
view(90,0)
title('Behind')
ylabel('y (mm)')
zlabel('z (mm)')

%-- [3] Head volume + channel in x-y plan
subplot(223)
checkreg_plot(vol, chanpos)
view(90,90)
title('Above')
xlabel('x (mm)')
ylabel('y (mm)')


% Slightly shift subplot to bottom
for i = 1:3
    subplot(2,2,i)
    posa = get(gca,'position');
    set(gca,'position',[posa(1) posa(2)-0.03 posa(3:4)])
end

titfig = {'MEG channels & head model obtained from segmented MRI:';infodata};
annotation(gcf,'textbox','String',titfig,'interpreter','none',...
    'FontSize',12,'fontname','AvantGarde','fontweight','bold',...
    'LineStyle','none','HorizontalAlignment','center',...
    'FitBoxToText','off','Position',[0.05 .88 0.9 0.12]);

if sav
    namfig = 'checkreg_chanvol';
    export_fig([savpath,filesep,namfig,'.jpeg'],'-m1.5')
    close
end

function checkreg_plot(vol, chanpos)

ft_plot_vol(vol, 'facecolor', [1 0.97 0.92], 'edgecolor', [0.46 1 0.73])
axis on
grid on
hold on
p = plot3(chanpos(:,1),chanpos(:,2), chanpos(:,3),'o');
set(p,'markersize',4,'markerfacecolor', [0.043 0.52 0.78])