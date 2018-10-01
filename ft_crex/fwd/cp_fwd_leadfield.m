function Sdb = cp_fwd_leadfield(Sdb)
% Prepare leadfield for subcortical and/or cortical sources depending on
% source.mat fieldnames (see cp_sources_*.m)
%
%   --> model.mat saved inside db_ft/PROJ/SUBJ/sources
%       containing headmodel, source grid with leadfield and associated brain atlas region 
%
% TO DO : CHECK FOR LEADFIELD normalize parameter
%
% See: % http://www.fieldtriptoolbox.org/tutorial/beamformer#compute_lead_field
%-- About leadfield normalization :
% "If you are not contrasting the activity of interest against another condition or baseline time-window, 
% then you may choose to normalize the lead field (cfg.normalize='yes'), which will help control
% against the power bias towards the center of the head."
%-- About channel selection:
% "Sensors MLP31 and MLO12 were removed from the data set. Thus it is 
% essential to remove these sensors as well when calculating the lead
% fields."
% --> MEG data preprocessing must be done before leadfield computation
% --> the resulting model.mat file is related to the MEG data that will be
% use for source analysis(cleanTrials.mat). It will be saved at the same
% location in meg/epoched/run_*/prep_dir (name of prep_dir depends on the
% filtering options).
%
%-CREx180530

Np = length(Sdb);

% Initialize waitbar
wb = waitbar(0, 'Leadfield computation...', 'name', 'Forward model');
wb_custcol(wb, [0 0.6 0.8]);

for i = 1 : Np
    psubj = Sdb(i);
    
    % Subject info
    sinfo = psubj.sinfo;
    
    rdir = psubj.meg.rundir;
    Nr = length(rdir);
    
	waitbar(i/Np, wb, ['Leadfield: ', sinfo]);
    for j = 1 : Nr
        % Check if computation already done according to MEG data
        % preprocessing (removing of bad channels)
        pmod = [psubj.meg.clean_dir{j}, filesep, 'fwd_model.mat'];
        % Already done
        if exist(pmod, 'file')
            Sdb(i).fwd.model_run{j} = pmod;
            continue
        end
        srun = rdir{j};
        waitbar(i/Np, wb, ['Leadfield: ', sinfo, '--', srun]);
        % Raw MEG data directory
        praw = psubj.meg.continuous.raw{j};
        pmeg = psubj.meg.clean_mat{j};
        if isempty(pmeg)
            warning('Preprocessed MEG data required for leadfield computation...')
            warning('Computation abort for subject %s\n', [sinfo, '--', srun]);
            continue
        end
        % Preproc MEG data directory
        prep = fileparts(pmeg);
        % Set forward model for both cortical and subcortical
        fwd_model = set_model(psubj, praw, prep);  

        % Save model data
        pmod = [prep, filesep, 'fwd_model.mat']; 
        save(pmod, 'fwd_model')

        % Update Sdb paths
        Sdb(i).fwd.model_run{j} = pmod;

        % Leadfield matrix figures
        pso = make_dir([prep, filesep, 'fwd_model']);
        % Cortical
        if ~isempty(fwd_model.cortical)
            leadfield_fig(fwd_model.cortical.grid, pso, 'cortical')
        end

        % Subcortical
        if ~isempty(fwd_model.subcortical)
            leadfield_fig(fwd_model.subcortical.grid, pso, 'subcortical')
        end

        %-- Figure of model geometry       
        cp_fwd_model_fig(fwd_model, pso)
    end
      
end
close(wb);
% Define fwd model for cortical and (or) subcortical sources
% model.cortical and model.subcortical both with field:
% - grid: source model with associated leadfield matrices
% - headmodel: volume of conduction
% - atlas: atlas information to identify region of each source
function model = set_model(psubj, praw, prep)
%-- Load data

% Conduction volume for figure of coregistration
pshell = psubj.fwd.shell;
volshell = loadvar(pshell);

% Sources (subcortical for now) -- TO DO : LOOP with cortical too
pso = psubj.fwd.sources;
sources = loadvar(pso);

% Grad
hdr = loadvar([praw, filesep, 'hdr_event.mat']);
Sgrad = hdr.grad;

% Get the final channel selection according to preproc.mat file associated
% with the preprpocessing cleanTrials.mat
preproc = loadvar([prep, filesep, 'preproc.mat']);
%%% Only for MEG here !!!
chansel = chan_sel(preproc.rm.sens);

%----
% Create lead field = forward solution
model = [];
model.subcortical = prepare_model(sources.subcortical, volshell, Sgrad, chansel);
model.cortical = prepare_model(sources.cortical, volshell, Sgrad, chansel);


% Prepare all that is needed for source estimation / beamforming: head
% model, source point grid and orientation, associated leadfield
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


% Keep source label
mdl = [];
mdl.grid = ldf_grid;
mdl.headmodel = shell;
mdl.atlas = grid.label;

% Figure showing the leadfield
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
    ylim([0.5 Nc+0.5])
    xlim([0.5 Ns+0.5])
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
