function winart = cmeg_artefact_detc(ftData, opt)
%--- Detect strong artefact based on Fieldtrip method
% http://www.fieldtriptoolbox.org/tutorial/automatic_artifact_rejection/
% ft_artifact_zvalue
% 
%-CREx181213

winart = [];

dopt = [];
dopt.rmchan = [];
% z-threshold for detection
dopt.zthr = 80;
% box car width for smoothing using ft_preproc_smooth (in second)
dopt.dt_smo = 0.200;
dopt.info = '';
dopt.zthr_low = 40;
% Duration for padding around artefact (left and right side) in second
dopt.dur_pad = 0.400;
dopt.force = 0;
% Maximum number of artefact (the threshold will be increase if the maximum
% number is exceeded
dopt.nart_max = 8;

if nargin < 2
    opt = dopt;
else
    opt = check_opt(opt, dopt);
end

% Detect artifact
Sart = prepare_detect(ftData, opt);
Sart = detc_artefact(Sart, opt.zthr);
if ~opt.force && isempty(Sart.ind)
    return;
end

% Not expected to remove more than opt.nart_max artefacts (default: 12)
% This can be due to too noisy dataset that would have to be excluded for
% analysis (not a valid run) --> during the bad sensor selection step
% Increase the threshold to have only 12 artefacts
if Sart.N > opt.nart_max
    Sart = adjust_thr(Sart, opt.nart_max);
end

% Message box
if ~opt.force   
    Na = Sart.N;
    if Na > 1
        ad = 's have';
    else
        ad = ' has';
    end
    uiwait(msgbox({['\fontsize{12}', num2str(Na) ' strong artefact', ad,' been detected'];...
        ['\fontsize{13}\bf ', prep_tex(opt.info)]}, 'Strong artefact', 'help',...
        struct('WindowStyle', 'non-modal', 'Interpreter', 'tex')));
else
    uiwait(msgbox({'\fontsize{12}Please adjust the threshold to detect strong artefacts';...
    ['\fontsize{13}\bf ', prep_tex(opt.info)]}, 'Strong artefact', 'help',...
    struct('WindowStyle', 'non-modal', 'Interpreter', 'tex')));
end

% Figure of artifact
artifact_fig(Sart);
% Add controls
add_buttons(Sart);

uiwait();

ia = ~isnan(get(findobj(gcf, 'tag', 'xart'), 'YData'));
if any(ia)
    iart = ind_art(ia);
    winart = ftData.time{1}(iart);
    winart = round(winart.*1e3)./1e3;
end
  
delete(gcf)

% Compute the z-sum across sensors
function Sart = prepare_detect(ftData, opt)

zthr = opt.zthr;
dt_smo = opt.dt_smo;
% Remove the bad channels first
time = ftData.time{1}; 
xd  = ftData.trial{1};  
label = ftData.label;

if ~isempty(opt.rmchan)
    isz = ismember(label, opt.rmchan);
    xd = xd(~isz, :);
    label = label(~isz);
end

[nc, nt] = size(xd);

% Remove the mean
xd = xd - repmat(mean(xd, 2), 1, nt);
% Envelope
xenv = abs(transpose(hilbert(transpose(xd))));

% Smooth envelope using ft function
Fs = fsample(time);
N = round(dt_smo * Fs);
if ~rem(N,2)
  % the kernel should have an odd number of samples
  N = N+1;
end
xenv = ft_preproc_smooth(xenv, N);
% compute the average and the standard deviation
% avg accross channels
menv = repmat(mean(xenv, 2), 1, nt);
sdenv = repmat(std(xenv, 0, 2), 1, nt);

% z-normalized envelope
zenv = (xenv - menv)./sdenv;

% accumulate the z-values over channels
zsum = sum(zenv)./sqrt(nc);

% Sort the channels depending on the abs amplitude at a zsum peak which is not
% at the border (peak due to the filtering)
% Excluding a 2s- border zone
int = time > 2 & time < time(end)-2;
tpk = time(int);
[~, izpk] = max(zsum(int));
ts = tpk(izpk);
[~, im] = sort(xenv(:, time==ts), 'descend');
zlab = label(im);

xdat = xd(im, :);
npad = round(Fs*opt.dur_pad);

Sart = [];
Sart.zsum = zsum;
Sart.time = time;
Sart.xdat = xdat;
Sart.xart = NaN(nc, nt);
Sart.chan = zlab;
Sart.thr = zthr;
Sart.npad = npad;
Sart.N = 0;
Sart.ichan = [1 nc];
Sart.ind = [];
Sart.thr_low = opt.zthr_low;

% Fieldtrip method but outside the ft_artifact_zvalue (faster)
function Sart = detc_artefact(Sart, zthr)
if ~isempty(Sart.ind) && Sart.thr==zthr
    return;
end
Sart.thr = zthr;
[nc, nt] = size(Sart.xart);
% threshold the accumulated z-values
art = Sart.zsum > zthr;

if ~any(art)
   Sart.ind = [];
   Sart.xart = NaN(nc, nt);
   Sart.N = 0;
   return;
end

% Extend at left and right depending on thr_low before padding
art = extend_art(Sart, art);

% pad the artifacts
npad = Sart.npad;
arti = find(diff([0 art])== 1)';
artf = find(diff([art 0])==-1)';
iart = [arti-npad artf+npad];
iart(iart < 0) = 1;
iart(iart > nt) = nt;

Na = length(iart(:, 1));
for i = 1 : Na
    int = iart(i, 1) : iart(i, 2);
    art(int) = 1;
end

iart = ind_art(art);
% Merge artefacts that are closed (< 2*npad)
[iart, art, Na] = fill_gap(iart, art, npad);

xn = NaN(nc, nt);
xn(:, art) = Sart.xdat(:, art);

Sart.ind = iart;
Sart.xart = xn;
Sart.N = Na;

function art = extend_art(Sart, art)
if Sart.thr <= Sart.thr_low
    return;
end
zthr = Sart.thr_low;
arti = find(diff([0 art])== 1)';
artf = find(diff([art 0])==-1)';
nt = length(art);
iart = [arti artf];
iart(iart < 0) = 1;
iart(iart > nt) = nt;
Na = length(iart(:, 1));
for i = 1 : Na
    ai = iart(i, 1);
    af = iart(i, 2);
    exti = find(Sart.zsum(ai:-1:1) <= zthr, 1, 'first');
    % Border cases
    if isempty(exti)
        exti = ai-1;
    end
    extf = find(Sart.zsum(af:end) <= zthr, 1, 'first');
    if isempty(extf)
        extf = length(Sart.zsum) - af;
    end
    int = ai-exti : af+extf;
    art(int) = 1;
end

% Merge artefacts that are closed (< 2*npad)
function [iart, art, Na] = fill_gap(iart, art, npad)
Na = length(iart(:, 1));
if Na > 1
    for i = 2 : Na
        if (iart(i, 1) - iart(i-1, 2)) <= 2*npad
            iart(i-1, 2) = iart(i, 2);
            iart(i, :) = [0 0];
        end
    end
    if any(iart(:,1)==0)
        iart = iart(iart(:, 1)~=0, :);
        Na = length(iart(:, 1));
        art = zeros(size(art));
        for i = 1 : Na
            int = iart(i, 1) : iart(i, 2);
            art(int) = 1;
        end
        art = logical(art);
    end
end
    
function Sart = adjust_thr(Sart, Nmax)
inc = 0.5;
tmax = 200;
while Sart.N >= Nmax && Sart.thr < tmax
    Sart = detc_artefact(Sart, Sart.thr + inc);
end

% Final artefacts including padding and no portion inside another
function iart = ind_art(isart)
arti = find(diff([0 isart])== 1)';
artf = find(diff([isart 0])==-1)';
iart = [arti artf];

% Figure showing dtected artefact + the z-sum profile for detection
function artifact_fig(Sart)

figure
set(gcf, 'units', 'centimeters', 'Position', [10 4 21 20]);
ax1 = axes(gcf, 'position', [0.08 0.60 0.87 0.29], 'Tag', 'ax_zsum');
ax2 = axes(gcf, 'position', [0.08 0.19 0.87 0.29], 'Tag', 'ax_data');

%--- Z-sum profile
plot_zsum(Sart);

%--- Data of one channel + artefacts
% First channel
set(gcf, 'UserData', Sart);
plot_data(Sart)

linkaxes([ax1, ax2], 'x')
h = zoom;
set(h, 'Enable', 'on', 'ButtonDownFilter', @notzoom)

annotation(gcf,'textbox','String','Validation of the strong artefact detection',...
    'FontSize',14, 'FontWeight', 'bold',...
    'LineStyle','none','HorizontalAlignment','Center',...
    'FitBoxToText','off','Position',[0.034 0.93 0.93 0.05],...
    'backgroundcolor',[1 1 1]);

% Plot new detection
function refresh(Sart)
set(gcf, 'UserData', Sart);
%--- Z-sum profile
plot_zsum(Sart);
%--- Data of one channel + artefacts
plot_data(Sart)

% Detcetion profile based on the z-normalized envelope signal averaged across
% channels
function plot_zsum(Sart)
ax = findobj(gcf, 'Tag', 'ax_zsum');
axes(ax);
hold off
time = Sart.time;
nt = length(time);
plot(time, Sart.zsum)
hold on
plot(time, ones(1, nt).*Sart.thr, 'color', [0.85 0 0], 'tag', 'thr_line');

xlim(time([1 end]))
grid on
set(gca, 'XMinorGrid', 'on', 'YMinorGrid', 'on')
lg = legend('zsum', 'threshold');
set(lg, 'fontsize', 13)
xlabel('Time (s)')
ylabel('z-sum')
set(gca, 'fontsize', 12)
title('Detection profile based on the normalized Z-average envelope of data across channels',...
    'fontsize', 13, 'FontWeight', 'normal')
set(ax, 'Tag', 'ax_zsum');

% Artefact highlight on the channel displaying max z value
function plot_data(Sart)
ax = findobj(gcf, 'Tag', 'ax_data');
axes(ax);
ichan = Sart.ichan;
ic = ichan(1);
hold off
time = Sart.time;
plot(time, Sart.xdat(ic, :));
hold on
plot(time, Sart.xart(ic, :), 'tag', 'xart')
yl = ylim;
xlim(time([1 end]))
set(gca, 'fontsize', 12)
ylim(yl)
put_rectangle(time(Sart.ind), [0.92 0.9 0.9])
xlabel('Time (s)')
ylabel('Amplitude')
title(['Detected artifacts highlighted on channel ', Sart.chan{ic}],...
    'fontsize', 13,'FontWeight', 'normal')
set(ax, 'Tag', 'ax_data');

function ino = notzoom(obj, ~)
% If the tag of the object is 'DoNotIgnore', then return true.
otag = obj.Tag;
if any(strcmp(otag, {'but_val'; 'slider'; 'ok_thr'; 'chan_prev'; 'chan_next'}))
   ino = true;
else
   ino = false;
end

% Launch the ft interactive detection
function new_detect(~, ~)
zthr = get(findobj(gcf, 'Tag', 'slider'), 'Value');
Sart = get(gcf, 'UserData');
Sart = detc_artefact(Sart, zthr);
refresh(Sart);

% Major tom
function add_buttons(Sart)

% Validate
uicontrol(gcf, 'Style', 'pushbutton',...
    'Units', 'normalized', 'Position', [0.61352 0.033121 0.22764 0.04586],...
    'BackgroundColor', [0.97 0.97 0.97], 'Fontsize', 13,...
    'ForegroundColor', [0 0.75 0],...
    'FontWeight', 'bold',...
    'String', 'Confirm detection',...
    'Tag', 'but_val',...
    'TooltipString', 'Reject the detected artefacts',...
    'Callback', @validate);

% Add panel to change channel
hchan = uipanel(gcf,'units', 'normalized',...
    'Position', [0.79663 0.48752 0.15283 0.062858],...
    'BackgroundColor', [1 1 1],...
    'title', 'Change channel',...
    'FontSize', 10);

uicontrol(hchan, 'Style', 'pushbutton',...
    'Units', 'normalized', 'Position', [0.19 0.195 0.24 0.60],...
    'BackgroundColor', [0.85 0.85 0.85], 'Fontsize', 9,...
    'FontWeight', 'bold',...
    'String', '<',...
    'Tag', 'chan_prev',...
    'TooltipString', 'Previous',...
    'Callback', {@chan_change, 'prev'});

uicontrol(hchan, 'Style', 'pushbutton',...
    'Units', 'normalized', 'Position', [0.575 0.195 0.24 0.60],...
    'BackgroundColor', [0.85 0.85 0.85], 'Fontsize', 9,...
    'FontWeight', 'bold',...
    'String', '>',...
    'Tag', 'chan_next',...
    'TooltipString', 'Next',...
    'Callback', {@chan_change, 'next'});

% Add a slider to modify threshol on the fly
hslid = uipanel(gcf,'units', 'normalized',...
    'Position', [0.27437 0.010191 0.25271 0.10446],...
    'BackgroundColor', [1 1 1],...
    'title', 'Change threshold',...
    'FontSize', 11);

uicontrol(hslid, 'Style', 'slider',...
    'Units', 'normalized',...
    'Position', [0.107 0.0812 0.396 0.377],...
    'Value', Sart.thr,...
    'Min', 0,...
    'Max', 400,...
    'SliderStep', [0.005 0.005],...
    'Tag', 'slider',...
    'Callback', @change_thr);

uicontrol(hslid, 'Style', 'text',...
    'Units', 'normalized', 'Position', [0.131 0.557 0.344 0.327],...
    'BackgroundColor', [1 1 1], 'ForegroundColor', [0.4 0 0], 'Fontsize', 12,...
    'FontWeight', 'bold',...
    'String', num2str(Sart.thr),...
    'Tag', 'thr_str');

uicontrol(hslid, 'Style', 'pushbutton',...
    'Units', 'normalized', 'Position', [0.689 0.295 0.208 0.409],...
    'BackgroundColor', [0.85 0.85 0.85], 'Fontsize', 12.5,...
    'FontWeight', 'bold',...
    'String', 'OK',...
    'Tag', 'ok_thr',...
    'TooltipString', 'Launch a new detection',...
    'Callback', {@new_detect});

set(gcf, 'CloseRequestFcn', []);

function chan_change(~, ~, sch)
% Get current number and max allowed number
Sart = get(gcf, 'UserData');
ichan = Sart.ichan;
if strcmp(sch, 'prev')
    ic = ichan(1) - 1;
    if ic < 1
        ic = ichan(2);
    end
else
    ic = ichan(1) + 1;
    if ic > ichan(2)
        ic = 1;
    end
end
ichan(1) = ic;
Sart.ichan = ichan;
set(gcf, 'UserData', Sart);
plot_data(Sart);

function change_thr(obj, ~)
thr = obj.Value;
set(findobj(gcf, 'Tag', 'thr_str'), 'String', num2str(thr));
hline = findobj(gcf, 'Tag', 'thr_line');
ydat = get(hline, 'YData');
set(hline, 'YData', ones(1, length(ydat))*thr, 'Color', [0.87 0.49 0]);

% Validate
function validate(~, ~) 
uiresume;

