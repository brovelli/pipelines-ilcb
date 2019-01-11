function cmeg_plot_data(ftData, opt)
% Figures showing continuous data recordings
% 6 channels are shown per figure with 3 temporal intervals: 
% - full (the entire duration) 
% - focus on a 10-s window from t0=100 s to t=110 s (to avoid some edge effects)
% - focus around the maximum absolute peak
% 

dopt = struct('savepath', [], 'info', [], 'colors', []);

opt = check_opt(opt, dopt);

% Check for directory to save figures
pfig = opt.savepath;
if isempty(pfig) || ~exist(pfig, 'dir')
    pfig = make_dir([pwd, filesep, 'datadisp'], 1);
end

Scol = opt.colors;
if ~isempty(Scol)
    colab = Scol.col;
    txtcol = repmat(mean(colab, 2) < 0.85, 1, 3);   
end
% Default colors for channel labels on plots
bcol = 'o';
tcol = 'w';

% Time
td = ftData.time{1};
td = td - td(1);

% Amplitude
xall = ftData.trial{1};

% Channels
label = ftData.label;

% Sort data by croissant channel number
numc = cellfun(@(x) str2double(regexprep(x, '\D*', '')), label);
[~, inds] = sort(numc);
label = label(inds);
xall = xall(inds, :);

% Channel number
Nc = length(label);
% Number of subplot per figure
nsub = 6;
vfig = 1 : nsub : Nc;
% Total number of figures
Nf = length(vfig);

% Vector to slightly shift subplot to the bottom
vbot = 0.008 : 0.005 : 0.008 + (nsub-1)*0.005;

xlz = [0 10] + 100; % Zoom sur le temps
tp = td(td>=xlz(1) & td<xlz(2));
tp = tp - tp(1);
xp = xall(:, td>=xlz(1) & td<xlz(2));
dtpeak = 4; % 4 s autour du pic max

wb = waitbar(0, 'Continuous data plots...', 'name', 'Cleanup figures');
wb_custcol(wb, [0 0.6 0.8]);
for i = 1 : Nf
    sn = num2str(i);
    snum = [repmat('0', 1, 2-length(sn)), sn];
    
    figure 
    set(gcf,'Visible','off','units','centimeters','position',[2 2 40 26])
    ch_ini = num2str(vfig(i));
    for j = 1 : nsub
        ichan = vfig(i) + j-1;
        ip = (j-1)*7;
        ysh = vbot(j);
        if ichan <= Nc
            chan = label{ichan};
            %------
            % Subplot des donnees MEG par canal - toute la duree
            subplot(nsub, 7, ip+1 : ip+3)
            plot(td(1:6:end), xall(ichan,1:6:end)) 
            format_ax([-0.05 -ysh 0.05 0]);
            
            % Add MEG channel name
            if ~isempty(Scol)
                bcol = colab(strcmp(Scol.chan, chan), :);
                tcol = txtcol(strcmp(Scol.chan, chan), :);
            end
            put_figtext(chan, 'nw', 12, tcol, bcol);
            % Titre si premiere ligne de subplot de la figure
            if j==1
                title('Continuous time series', 'fontweight','bold', 'interpreter', 'none')
            end
            
            %------
            % Subplot des donnees MEG par canal - 10 premieres secondes
            subplot(nsub,7, ip+4 : ip+5)
            plot(tp, xp(ichan, :))
            format_ax([0.005 -ysh 0.005 0]);           
            put_figtext(chan, 'nw', 12, tcol, bcol);
            % Titre si premiere ligne de subplot de la figure
            if j==1
                title(['10 s from t0=',num2str(xlz(1)),' s'],'fontweight','bold')
            end 
            
            %------
            % Subplot zoom autour du max d'amplitude
            subplot(nsub,7,ip+6:ip+7)
            plot_peak(td, xall(ichan, :), dtpeak)
            format_ax([0.022 -ysh 0 0]);   
            if j==1
                title('Around maximum amplitude','fontweight','bold')
            end           
        end
    end
    if ichan >= Nc
        iend = Nc;
    else
        iend = ichan;
    end
    ch_end = num2str(iend);
    %------
    % Titre general de la figure contenant le chemin d'acces aux donnees 
    tit={['Data display -  [', snum,']'] ; opt.info};
    
    annotation(gcf,'textbox','String',tit,'interpreter','none',...
        'FontSize',13,'fontname','AvantGarde',...
        'LineStyle','none','HorizontalAlignment','center',...
        'FitBoxToText','off','Position',[0.1 0.88 0.9 0.12]);

    namfig = ['datadisp_', snum,'_chan_', ch_ini,'_to_', ch_end];
    export_fig([pfig,filesep,namfig,'.png'], '-m1.2', '-nocrop')
    close
    waitbar(i/Nf);
end

close(wb);
fprintf('\n\nLook at figures in ---\n----> %s\n', pfig)

% Format axis 
function format_ax(shift)
    axis tight
    ylabel('Magnetic field (T)')
    xlabel('Time (s)')
    % Move subplot down
    set_pos(shift)
    
function set_pos(shift)
pos = get(gca,'position');
set(gca,'position', pos + shift)    

function plot_peak(td, xx, dtpeak)
% Avoid peak at edge (search for peak outside data border of 4 s) 
dtbor = 4;
tc = td(td>td(1)+dtbor & td<td(end)-dtbor); 
xc = xx(td>td(1)+dtbor & td<td(end)-dtbor);

if any(abs(xc) > 0)
    [~, imax] = max(abs(xc));  
else
    imax = round(length(xc)./2);
end
tbor = [tc(imax)-dtpeak tc(imax)+dtpeak];
xcp = xc(tc>=tbor(1) & tc<=tbor(2));
tcp = tc(tc>=tbor(1) & tc<=tbor(2));
plot(tcp, xcp);
