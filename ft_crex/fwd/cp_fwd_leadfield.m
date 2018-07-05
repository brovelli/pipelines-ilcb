function Sdb = cp_fwd_leadfield(Sdb)
% Prepare leadfield for subcortical and/or cortical sources depending on
% source.mat fieldnames (see cp_sources_*.m)
%
%   --> model.mat saved inside db_ft/PROJ/SUBJ/sources
%       containing headmodel, source grid with leadfield and associated brain atlas region 
%
% TO DO : CHECK FOR LEADFIELD normalize parameter
%
%-CREx180530

Np = length(Sdb);


for i = 1 : Np
    psubj = Sdb(i);
    % Check if already done

    if isfield(psubj, 'model') && exist(psubj.model, 'file')
        % Already done
        continue
    end
    
    % Set forward model for both cortical and subcortical
    model = set_model(Sdb(i));  

    % Save model data
    pfwd = fileparts(Sdb(i).shell);
    pmod = [pfwd, filesep, 'model.mat']; 
    save(pmod, 'model')

    % Update Sdb paths
    Sdb(i).model = pmod;

    % Leadfield matrix figures
    % Cortical
    if ~isempty(model.cortical)
        leadfield_fig(model.cortical.grid, pfwd, 'cortical')

    end

    % Subcortical
    if ~isempty(model.subcortical)
        leadfield_fig(model.subcortical.grid, pfwd, 'subcortical')
    end

    %-- Figure of model geometry       
    cp_fwd_model_fig(model, pfwd)
      
end

% Define fwd model for cortical and (or) subcortical sources
% model.cortical and model.subcortical both with field:
% - grid: source model with associated leadfield matrices
% - headmodel: volume of conduction
% - atlas: atlas information to identify region of each source
function model = set_model(psubj)
%-- Load data

% Conduction volume for figure of coregistration
pshell = psubj.shell;
volshell = loadvar(pshell);

% Sources (subcortical for now) -- TO DO : LOOP with cortical too
pso = psubj.sources;
sources = loadvar(pso);


% Path to raw MEG data %%%% TO DO : to preprocessed / epoched/run_concat
pmeg = psubj.meg.continuous.raw{1};

%- Looking for the raw MEG data in pmeg directory
draw = filepath_raw(pmeg);

% Get sensor position
Sgrad = ft_read_sens(draw);

%----
% Create lead field = forward solution

% Only select MEG label 
% Sgrad = Strials.(fcond{1}).grad;
% chanlab = Strials.(fcond{1}).label'; 
% % ftTuto : "essential to indicate removing sensors when calculating the lead fields"
% Sgrad = ft_convert_units(Sgrad,'mm');
%%% TO DO - chansel from preprocessing MEG data !!!
chansel = ft_channelselection('meg', Sgrad.label);

model = [];
model.subcortical = prepare_model(sources.subcortical, volshell, Sgrad, chansel);
model.cortical = prepare_model(sources.cortical, volshell, Sgrad, chansel);



function mdl = prepare_model(grid, shell, grad, chan)
if isempty(grid)
    mdl = [];
    return
end

% Prepare source model
% Add 'inside' field if necessary (if customised grid)
if ~isfield(grid, 'inside')
    cfg = [];
    cfg.grid = grid;
    cfg.headmodel = shell;
    grid = ft_prepare_sourcemodel(cfg);
end

if isfield(grid, 'mom')
    grid.mom = grid.mom';
end

cfg = [];
cfg.grad = grad;
cfg.headmodel = shell;
cfg.grid = grid;
cfg.reducerank = 2;
cfg.channel = chan;  
ldf_grid = ft_prepare_leadfield(cfg);

% Figure showing the leadfield

% cfg.normalize = 'yes';
% If you are not contrasting the activity of interest against another condition or baseline time-window, 
% then you may choose to normalize the lead field (cfg.normalize='yes'), which will help control
% against the power bias towards the center of the head. 

% Keep source label
mdl = [];
mdl.grid = ldf_grid;
mdl.headmodel = shell;
mdl.atlas = grid.label;

function leadfield_fig(ldf_grid, pfig, styp)

Nc = length(ldf_grid.label);
ins = ldf_grid.inside;
Ns = sum(ins);

ldfs = ldf_grid.leadfield{find(ins==1, 1, 'first')};
Nd = length(ldfs(1, :));

ldf = cell2mat(ldf_grid.leadfield);

ldf = reshape(ldf, Nc, Nd, Ns);

if Nd==1
    sdir = {'mom'};
else
    sdir = {'x', 'y', 'z'};
end

figure

set(gcf, 'visible', 'off', 'units','centimeters','position', [5 5 30 15])

for i = 1 : Nd
    
    subplot(1, 3, i)
    
    imagesc(squeeze(ldf(:, i, :)))

    xlabel('Dipole number','fontsize',12)
    ylabel('Channel number','fontsize',12)
    
    set(gca,'ydir','normal', 'fontsize',11, 'fontweight', 'bold')

    title({'Leadfield matrix';['in ', sdir{i}, '-direction']},...
        'fontsize',12, 'fontweight','normal')
end

colormap(colormap_blue2red)
caxis([min(ldf(:)) max(ldf(:))])

for i = 1 : Nd
    subplot(1, 3, i)
    pos = get(gca, 'position');
    set(gca, 'position', [pos(1)-0.006 pos(2:end)])
end
g = colorbar;
if Nd >= 3
    set(g, 'position', [0.9344 0.1664 0.0122 0.3201])
else
    set(g, 'position', [0.369 0.1664 0.0122 0.3201])
end

export_fig([pfig, filesep, 'leadfield_', styp,'.png'], '-m2', '-c[NaN,NaN,NaN,NaN]', '-p0.01')
close
