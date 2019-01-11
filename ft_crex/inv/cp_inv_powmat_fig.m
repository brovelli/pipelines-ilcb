function cp_inv_powmat_fig(Smean, pfig)
% Figures showing power accross trials for each ROI
%-CREx180920
if nargin < 2 || isempty(pfig)
    pfig = make_dir('powmat_fig', 1);
end

% Info for figure title
Sinf = Smean.info;
fqnam = Sinf.param.name;
cond = Sinf.cond_out;
info = Sinf.subj;

pow = Smean.pow;
time = Smean.time;

% Number of trials
ntr = length(pow(1, :, 1));
% Number of ROI
Na = length(Smean.label);

% Colormap limits based on median +/- a factor*standard deviation
md = nanmedian(pow(:));
sd = nanstd(pow(:));
mpow = [md-2*sd md+4*sd];

% Loop over ROI
for j = 1 : Na
    figure('visible', 'off')
    
    imagesc(time, 1:ntr, squeeze(Smean.pow(j, :, :)));
    colormap(colormap_spectral)
    xlim(time([1 end]))
    
    set(gca, 'fontsize', 12)
    caxis(mpow) 
    
    % Colorbar
    cb = colorbar;
    set(cb, 'position', [0.91726 0.11667 0.020238 0.32143])
    cb.Label.String = 'z-score';
    cb.Label.Position = [1.9152 0.5 0];
    cb.Label.FontWeight = 'bold';
    cb.Label.FontSize = 11;

    tit = {['[', fqnam, '] ', info, ' - Cond: ', cond]; Smean.labfull{j}};
    title(tit , 'fontsize', 14, 'interpreter', 'none')
    ylabel('Trials')
    xlabel('Time (s)')
    hold on
    line([0 0], [0 ntr], 'linestyle', ':', 'linewidth', 0.8, 'color', 'k');
    export_fig([pfig, filesep, 'pow_', Smean.label{j}, '_', fqnam,'_', cond,'.png'], '-m1.2')
    close
end  