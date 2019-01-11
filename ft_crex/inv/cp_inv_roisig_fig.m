function cp_inv_roisig_fig(time, csig, opt)
%
%-CREx-190110
dopt = struct('name', [],...
            'ci_sig', [],...
            'meshes', [],...
            'dip_pos', [],...
            'dip_ori', [],...
            'dip_col', [],...
            'ylim_abs', [],...
            'side', 'L',...
            'title', '',...
            'sig_title', 'Mean power',...
            'info', '',...
            'savepath', pwd,...
            'visible', 'on');
            
opt = check_opt(opt, dopt);

if ~iscell(csig)
    csig = {csig};
end

ci_sig = opt.ci_sig;
if ~isempty(ci_sig) && ~iscell(ci_sig)
    ci_sig = {ci_sig};
end
    
if isempty(opt.meshes)
    warning('No meshes in input parameter (opt structure)')
end

if ~isempty(opt.info)
    sinfo = ['_', opt.info];
else
    sinfo = '';
end

%- Directory to save figure
pfig = opt.savepath;

%- Define dipole spheres
dip_pos = opt.dip_pos;
Sbals = dipole_spheres(dip_pos, [], 2, opt.dip_col);
% Orientation
dip_ori = opt.dip_ori;

meshes = opt.meshes;

%- Number of signal to superimpose
Ns = length(csig);
%- Color of each signal
colp = color_group(Ns);

% Defin a relative ylim
% yrel = [cellfun(@min, csig, 'UniformOutput', 1) cellfun(@max, csig, 'UniformOutput', 1)];
% yrel = [min(yrel(:, 1)) max(yrel(:, 2))];
% dy = diff(yrel)./100;
% yrel = [yrel(1)-dy yrel(2)+dy];

%- Absolute lim
yabs = opt.ylim_abs;
if ~isempty(yabs)
    dy = diff(yabs)./100;
    yabs = [yabs(1)-dy yabs(2)+dy];
    pfiga = make_dir([pfig, filesep, 'scale_absolute']);
end
pfigr = make_dir([pfig, filesep, 'scale_relative']);

%---- Figure
figure('visible', opt.visible, 'units', 'centimeters',...
    'position', [3 6 24 15], 'color', [0 0 0])

%- Mesh view 1
axg_1 = axes(gcf, 'position', [0.02 0.5 0.34 0.38], 'color', [0 0 0]);
%- Mesh view 2
axg_2 = axes(gcf, 'position', [0.02 0.06  0.34 0.38]);
%- Signal plot
ax_sig = axes(gcf, 'position', [0.45 0.25 0.5 0.42]);

%- Plot MESH 1
axes(axg_1);
mopt = [];
mopt.newfig = 0;
mopt.dispnorm = 0;
[hp, hleg] = plot_meshes(meshes, mopt);
delete(hleg);

% Change mesh color and alpha
set(hp, 'edgecolor', 'none');
set(hp(1), 'facealpha', 0.13)
if length(hp) > 1
    set(hp(2), 'facealpha', 0.12)
end
hold on

% Add dipoles
arrayfun(@patch, Sbals)

% Add norms
if ~isempty(dip_ori)
    hq = quiver3c(dip_pos, dip_ori);
    set(hq, 'AutoScaleFactor', 1.3, 'linewidth', 0.2, 'color', [0.7 0.2 0], 'PickableParts', 'none');
else
    hq = [];
end

% View-1
view(0, 90)
camlight('right','infinite')
axis tight equal off;
ao = findobj(gca, 'type', 'patch');

axes(axg_2)
hold on
copyobj([ao; hq], axg_2)

%- View-2
if strcmpi(opt.side, 'L')
    view(180, 0);
else
    view(0, 0);
end

camlight('right','infinite')
axis tight equal off;

%- Plot signals
axes(ax_sig)
set(gca, 'xcolor', [1 1 1], 'ycolor', [1 1 1])
hold on
hsig = zeros(Ns, 1);
for i = 1 : Ns
    hsig(i) = plot(time,  csig{i}, 'color', colp(i, :), 'linewidth', 1.2);
    if ~isempty(ci_sig)
        [hb, pci] = boundedline(time, csig{i}', ci_sig{i}');
        set(pci, 'facecolor', colp(i, :), 'facealpha', 0.55);
        set(hb, 'linewidth', 0.8, 'color', colp(i,:));
    end
end
xlim(time([1 end]))
% ylim(yrel)

set(gca, 'color', [0 0 0]);
grid on
box on

set(gca, 'GridColor', [0.92 0.92 0.92], 'GridAlpha', 0.4)
title(opt.sig_title, 'fontsize', 12, 'color', [1 1 1])
xlabel('Time (s)', 'color', [1 1 1])
ylabel('Normalized power (z-score)', 'color', [1 1 1])
set(gca, 'fontsize', 11)

if ~isempty(opt.name)
    hleg = legend(hsig, opt.name, 'location', 'northwest');
    set(hleg, 'position', [0.84304 0.077359 0.11355 0.082011],...
        'fontsize', 12, 'color', 'none',...
        'TextColor', [1 1 1], 'EdgeColor', 'none');
end

annotation(gcf, 'textbox', [0.101 0.91998 0.79926 0.061736], 'color', [1 1 1],...
    'fontweight', 'bold', 'fontsize', 14, 'String', opt.title,...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle',...
    'interpreter', 'none')
yl = ylim;
hl = line([0 0], yl, 'linestyle', '--', 'color', [1 1 1]);
ylim(yl);
export_fig([pfigr, filesep, 'roisig', sinfo, '_rel.png'], '-m1.5')

if ~isempty(yabs)
    delete(hl)
    line([0 0], yabs, 'linestyle', '--', 'color', [1 1 1]);
    ylim(yabs)
    export_fig([pfiga, filesep, 'roisig', sinfo, '_abs.png'], '-m1.5')
end
close    