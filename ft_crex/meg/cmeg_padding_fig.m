function cmeg_padding_fig(padData, iniData, opt)
% Make figure of continuous data to see the padding artefact effects
% Subplot of the padding data vs the original data, for the 3 channels that
% displayed the most important amplitude before the artefact correction
%
%-CREx180726

dopt = struct('savepath', [], 'info', []);
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
if ~isfield(padData, 'artpad') || isempty(padData.artpad)
    return;
end

Na = length(padData.artpad(:, 1));

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
xdat = iniData.trial{1};
xpad = padData.trial{1};

% Plot the data at the first 3 channels displaying the max amplitude before
% padding process
[~, ind] = sort(max(abs(xdat),[],2), 'descend');

xl = time([1 end]);

for i = 1:3
    k = ind(i);
    chan = iniData.label{k};
    tit{1} = ['Padding-artefact result for ', chan];
    
    figure
    set(gcf,'units','centimeter','position',[1 1 27 17])

    %- Padding data
    subplot(211) 
    plot(time, xpad(k, :))
    title(tit, 'fontsize', 14, 'interpreter', 'none')   
    format_ax(xl)

    %- Original data    
    subplot(212)
    plot(time, xdat(k, :))
    format_ax(xl)

    put_figtext('B) Original data','nw', 12, 'w', 'o');
    
    pos = get(gca,'position');
    set(gca,'position',[pos(1) pos(2)+0.05 pos(3:4)])

    subplot(211), put_figtext('A) Padded artefacts','nw',12, 'w', 'o');
    
    padnam = ['padart_', srma, '_', chan];
    
    export_fig([pfig, filesep, padnam, '.png'], '-m1.5')
    close
end


function format_ax(xl)

xlabel('Time (s)','fontsize',14)
ylabel('Amplitude','fontsize',14)
set(gca,'fontsize',14)
xlim(xl)