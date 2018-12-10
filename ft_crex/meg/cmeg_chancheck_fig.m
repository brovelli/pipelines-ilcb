function Scol = cmeg_chancheck_fig(ftData, savpath, info)
% Display the layout of MEG channels with an indication about the mean
% value of the Hilbert envelop of the data at each channel
% The value is normalized by the minimum enveloppe value found across all 
% the channels. 
% A value of 4 for channel A42 indicate that the mean amplitude of Hilbert
% envelop of the data recording by this channel is 4 order higher than the
% minimum mean value found in the whole data set.

if nargin<3
    info = '';
end
if nargin<2
    savpath = pwd;
end

if ~isfield(ftData, 'grad')
    disp('No grad field in data (probably SEEG data)')
    disp('Abort the diplay of the layout with mean envelop amplitudes')
    Scol = [];
    return;
end
chan = ftData.label;
Nc = length(chan);

% Define if 4D or neuromag data
if length(chan{1})> 3 && strcmp(chan{1}(1:3), 'MEG')
    flay = 'neuromag306all.lay';
else
    flay = '4D248_helmet.mat';
end

% Remove mean
xd = ftData.trial{1};
xd = xd - repmat(mean(xd,2),1,length(xd(1,:)));

% Only channels with values ~=0
isok = sum(xd, 2)~=0;

th_1 = hilbert(xd);               
th_2 = sqrt(xd.^2 + th_1.*conj(th_1));
xmh = mean(th_2, 2);

% Define colormap
cmap = colormap_blue2red;
val = linspace(min(xmh(isok)), max(xmh), length(cmap(:,1)));

% Get the x,y sensor position on layout
cfg = [];
cfg.grad = ftData.grad;
cfg.layout = flay;
lay = ft_prepare_layout(cfg);
lay.label = lay.label(1:end-2);
lay.pos = lay.pos(1:end-2,:);

xy = lay.pos;

% Channels that have been removed from analysis because of important noise...
laylab = lay.label;
[isc, indc] = ismember(laylab, chan);
if Nc < length(laylab) 
    xy = xy(isc, :);
    laylab = laylab(isc);
end
indc = indc(indc~=0);
xmh = xmh(indc);
isok = isok(indc);

chancol = cell2mat(arrayfun(@(x) cmap(find(val <= x, 1, 'last'), :), xmh, 'UniformOutput', 0));

if any(~isok)
    chancolz = zeros(Nc, 3);
    chancolz(isok, :) = chancol;
    chancolz(~isok, :) = repmat([0.4 0.4 0.4], sum(~isok), 1);
    chancol = chancolz;
end

% Keep color data for time series disp figures
Scol = [];
Scol.chan = laylab;
Scol.col = chancol;
Scol.val = xmh;
Scol.pos = xy;

% Plot
hb_layplot(lay, Scol, cmap);

tit = {info; 'Mean value of signal envelop (from hilbert transform)'};

annotation(gcf,'textbox','String',tit,'interpreter','none',...
    'FontSize',11,'fontname','AvantGarde',...
    'LineStyle','none','HorizontalAlignment','Left',...
    'FitBoxToText','off','Position',[0.0033 0.9427 0.9489 0.0525],...
    'backgroundcolor',[1 1 1]);

pfig = [savpath, filesep, 'chancheck_', num2str(Nc), 'chan.png'];
export_fig(pfig, '-m1.5')
close
fprintf('\nFigure showing mean values saved as :\n %s\n', pfig)

function hb_layplot(lay, Scol, cmap)
xy = Scol.pos;
txtcol = repmat(mean(Scol.col, 2) < 0.85, 1, 3);
figure
set(gcf, 'visible', 'off', 'units','centimeters','position',[3 1 24 24])    
set(gca,'position',[0 0.005 1 0.95])
ft_plot_lay(lay, 'point', true, 'pointsymbol','.','pointsize',1,'pointcolor',[1 1 1],...
    'box', false, 'label', true, 'labelsize',8,'mask', false, 'outline', true);

hold on

Nc = length(Scol.chan);
for k = 1 : Nc
    tx = text(xy(k, 1), xy(k, 2), Scol.chan{k});
    set(tx,'fontsize',9,'fontweight','bold','backgroundcolor', Scol.col(k, :),'color', txtcol(k, :))
end
colormap(cmap);
cb = colorbar;

set(cb,'position', [0.91 0.024 0.029 0.186])
yl = get(cb,'ylim');
set(cb,'ytick',[yl(1) (yl(2)-yl(1))/2 yl(2)])
val = Scol.val;
mi = min(val);
ma = max(val);
set(cb,'yticklabel',{'1',num2str((ma-mi)./(2*mi),'%3.1f'),...
    num2str(ma./mi,'%3.1f')})  
set(get(cb,'ylabel'),'string','Normalized value')

