function cmeg_artefact_fig(padData, iniData, opt)
% Subplot showing the padding and the original data, for the 3 channels that
% were displaying the most important amplitude before the artefact correction
%
%-CREx180726

dopt = struct('savepath', [], 'info', [], 'rmchan', []);
if nargin < 3 
    opt = dopt;
else
    opt = check_opt(opt, dopt);
end
if isfield(opt, 'figpath') && ~isempty(opt.figpath)
    opt.savepath = opt.figpath;
elseif isempty(opt.savepath)
    opt.savepath = pwd;
end

% Number of artefact to remove
if ~isfield(padData.hdr, 'artpad') || isempty(padData.hdr.artpad)
    return;
end
iwin = padData.hdr.artpad;

% Exclude bad channels
if ~isempty(opt.rmchan)
	cfg = [];
	cfg.channel = chan_sel(opt.rmchan);
	padData = ft_preprocessing(cfg, padData);
	iniData = ft_preprocessing(cfg, iniData);
end

Na = length(iwin(:, 1));

% Figure directory
srma = ['padart_', num2str(Na), 'ra'];
if isempty(opt.savepath)    
    pfig = make_dir([pwd, filesep, srma], 1);
else
    pfig = make_dir(opt.savepath);
end

tit = cell(2,1);
tit{2} = opt.info;

time = iniData.time{1};
labi = iniData.label;

xini_all = iniData.trial{1};
labp = padData.label;
xpad_all = padData.trial{1};

% Plot the data at the first 3 channels displaying the max amplitude before
% padding process
[~, ind] = sort(max(abs(xini_all),[],2), 'descend');
labsel = labp(ind(1:3));
Nc = length(labsel);
nt = length(time);
for i = 1 : Nc
    
    chan = labsel{i};
    
    xini = xini_all(strcmp(labi, chan), :);
    xpad = xpad_all(strcmp(labp, chan), :);
    % Highlight artifacts portions
    xred = NaN(1, nt);
    % Corrected portions
    xcor = NaN(1, nt);
    
    ipad = (xini - xpad) ~= 0;
    xred(ipad) = xini(ipad);
    xcor(ipad) = xpad(ipad);
    
    tit{1} = ['Artefact rejection for ', chan];
    
    figure
    set(gcf, 'visible', 'off', 'units', 'centimeter', 'position', [1 1 27 17])

    %- Padding data
    subplot(211) 
    plot(time, xpad)
    hold on
    plot(time, xcor)
    title(tit, 'fontsize', 14, 'interpreter', 'none')   
    format_ax(time)
    put_rectangle(time(iwin), [0.92 0.9 0.9])
    
    %- Original data    
    subplot(212)
    plot(time, xini)
    hold on
    plot(time, xred)
    format_ax(time)
    put_rectangle(time(iwin), [0.92 0.9 0.9])
    
    put_figtext('B) Original','nw', 12, 'w', 'o');
    
    pos = get(gca,'position');
    set(gca,'position',[pos(1) pos(2)+0.05 pos(3:4)])

    subplot(211), put_figtext('A) Artefact padding','nw',12, 'w', 'o');
    
    padnam = ['rmart_', srma, '_', chan];
    
    export_fig([pfig, filesep, padnam, '.png'], '-m1')
    close
end

function format_ax(time)
yl = ylim;
xlabel('Time (s)','fontsize',14)
ylabel('Amplitude','fontsize',14)
set(gca,'fontsize',14)
xlim(time([1 end]))
ylim(yl);