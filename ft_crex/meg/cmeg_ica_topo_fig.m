function cmeg_ica_topo_fig(compData, opt)
% Figures to show ICA components
%
% For each component : 
%       - time signal for the whole recording duration
%       - a detail of the time signal in a shorter time windows
%       (see opt.xlim_detail as input parameter)
%       - a topographic representation (spatial interpolation on the MEG layout)
% 
% If data are not MEG data (i.e. SEEG), the topographic representation is not performed.
%
%-------------------
% Input parameters
%
% -- compData : the Fieldtrip structure as output by ft_componentanalysis
% -- opt : structure of parameters with fields :
%
%       opt.savepath : directory path for saving figures [default: pwd/ica_fig]
%
%       opt.info : additionnal information to be added to the figures titles
%       [default: ''] (ex. the path to the component data set)
%
%       opt.xlim_detail : time window to display component signal detail (in second)
%       [default : [100 110] s] - to deactivate this option, set opt.xlim_detail
%       to 'no' or to empty ([])
%
%       opt.mode :  - 'auto' (figure are invisible and automatically print in
%                   savepath directory) [default]
%                   - 'check' : only figure of time signal + component are shown
%                   for a manual check of the signal (figure are not printed)
%
%--------------------
%
% This function use the fieldtrip function ft_topoplotIC for the topographic
% plots and function from ft_crex (including export_fig package)
%
%-CREx-180424

% Default options
dopt = struct('info', '',...
    'savepath', '',...
    'xlim_detail', [100 110],...
    'mode', 'auto');

% Check for inputs
if nargin < 2
    opt = dopt;
else
    opt = check_opt(opt, dopt);
end

% Define layout according to label
label = compData.topolabel;

% Define if 4D or neuromag data
[flay, sun] = cmeg_det_lay(label);

% Figure directory
if isempty(opt.savepath) || ~exist(opt.savepath, 'dir')
    pfig = make_dir([pwd, filesep, 'ica_fig'], 1);
else
    pfig = opt.savepath;
end

% Mode
if strcmp(opt.mode, 'check')
    check = 1;
else
    check = 0;
end

% Additionnal figure info 
sinfo = opt.info;

% Number of time signal per figure for combined compo + time signal plots
Nsub = 5;

% Check for detailed time window
time = compData.time{1};
if isempty(opt.xlim_detail) || ischar(opt.xlim_detail)
    xl = [];
else
    xl = opt.xlim_detail;
    if length(xl)==1
        xl = [xl xl+10];
    end
    ti = time(1); % == 0
    tf = time(end);
    if tf <= 15 && xl(2)>=15
        xl = [];
    else
        if xl(1) > tf || xl(1) < ti || xl(2) > tf || xl(2) < ti
            xl = [5 15];                
        end
    end
end

if ~check  
    toposig_plot(compData, flay, xl, Nsub, sinfo, sun, pfig)    

else
    toposig_check(compData, flay, Nsub, sinfo, sun)
end

% Topographic + time signal (Nsub components per figure)
function toposig_plot(compData, flay, xzoom, Nsub, sinfo, sun, pfig, check)
if nargin < 8
    check = 0;
end

if check==1
    svis = 'on';
else
    svis = 'off';
end
if ~isempty(sinfo)
    sinfo = ['data: ', sinfo];
end
 
sylab = ['Amplitude (', sun, ')'];

% ICA method
meth = compData.cfg.method;

% Number of components
Nc = length(compData.label);

subc = 1 : Nsub : Nc;
Nf = length(subc);

time = compData.time{1};
xcomp = compData.trial{1};

xl_full = time([1 end]);
xl_zoom = xzoom;

for i = 1 : Nf 
    
    ideb = subc(i);
    
    sn = num2str(i);
    if length(sn)==1
        sn = ['0', sn]; %#ok
    end
    
    figure 
    set(gcf, 'visible', svis, 'units', 'centimeters', 'position', [7 3 34 25])
    
    hax = zeros(1, Nsub);
    for j = 1 : Nsub
        k = ideb + j-1;

        if k <= Nc   
            if flay
                subplot(Nsub, 4, 4*j-3)
                cfg = [];
                cfg.component = k;
                cfg.layout    = flay;
                cfg.comment   = 'no'; 
                ft_topoplotIC(cfg, compData);
                pos = get(gca,'position');
                set(gca,'position', [0.05 pos(2:4)])
            end         
            % Time signal
            hax(j) = subplot(Nsub, 4, 4*j-2 : 4*j);
            plot(time, xcomp(k,:));
            xlim(xl_full)
            pos = get(gca,'position');
            set(gca,'position', [pos(1)-0.08 pos(2) pos(3)+0.1 pos(4)])
            ylabel(sylab, 'fontsize', 12)
            set(gca, 'fontsize', 12)
            put_figtext(k, 'nw', 13, [1 1 1], [0.8 0.3 0.1]);
        end
    end
    xlabel('Time (s)', 'fontsize', 12)
    if k >= Nc
        iend = Nc;
    else
        iend = k;
    end
    
    tit = {['ICA (', meth, ' method) - Components plot [', sn,']'] ; sinfo};
    
    annotation(gcf,'textbox','String', tit, 'interpreter','none',...
        'FontSize',13,'fontname','AvantGarde',...
        'LineStyle','none','HorizontalAlignment','center',...
        'FitBoxToText','off','Position',[0.1 0.88 0.9 0.12]);

    if ~check
        % Save figure    
        namsav = ['composig_',sn,'_comp_', num2str(ideb), '_to_', num2str(iend)];
        export_fig([pfig, filesep, namsav,'.png'],'-m1.5', '-c[20,20,20, NaN]')
        
        if ~isempty(xl_zoom)
            for j = 1 : Nsub
                if hax(j)~=0
                    set(hax(j), 'xlim', xl_zoom);
                end
            end
            export_fig([pfig, filesep, namsav,'_detail.png'],'-m1.5', '-c[20,20,20, NaN]')  
        end
    else
        zoom on
        linkaxes(hax(hax~=0), 'x')
        goon = input('See next components (1) or stop (0) ? ');
        if ~goon
            close
            break;
        end
    end    
    close
end

function toposig_check(compData, flay, Nsub, sinfo, sun)
toposig_plot(compData, flay, [], Nsub, sinfo, sun, [], 1)