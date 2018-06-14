function Sdb = cp_fwd_leadfield(Sdb)
% Prepare leadfield for subcortical and/or cortical sources depending on
% source.mat fieldnames (see cp_sources_*.m)
%
%   --> fwdmodel.mat saved inside db_ft/PROJ/SUBJ/sources
%       containing head model, source and leadfield for subsortical regions
%
% TO DO : ADD CORTICAL SOURCES (only subcortical here)
%-CREx180530

Np = length(Sdb);


for i = 1 : Np
    psubj = Sdb(i);
    % Check if already done
    adone = 0;
    if isfield(psubj, 'model')
        pmo = psubj.model;
        model = loadvar(pmo);
        % Only subcortical for now !!
        if isfield(model, 'subcortical') && ~isempty(model, 'subcortical')
            adone = 1;
        end
    end
    
    if ~adone
        
        modl = set_model(Sdb(i));  %--TO DO - include cortical sources !!
        
        model.subcortical = modl.subcortical;
    
        % Save model data
        pfwd = fileparts(Sdb(i).shell);
        pmod = [pfwd, filesep, 'model.mat']; 
        save(pmod, 'model')
        Sdb(i).model = pmod;

        %-- Figure of model geometry       
        model_geom_fig(modl.subcortical, pfwd)
    end
        
end

function model = set_model(psubj)
%-- Load data

% Conduction volume for figure of coregistration
pshell = psubj.shell;
shell = loadvar(pshell);

% Sources (subcortical for now) -- TO DO : LOOP with cortical too
pso = psubj.sources;
sources = loadvar(pso);
subso = sources.subcortical;

% Path to raw MEG data
pmeg = psubj.meg;

%- Looking for the raw MEG data in pmeg directory
draw = filepath_raw(pmeg);

% Get sensor position
Sgrad = ft_read_sens(draw);

% Prepare source model
% Add 'inside' field if necessary (if customised grid)
if ~isfield(subso, 'inside')
    cfg = [];
    cfg.grid = subso;
    cfg.headmodel = shell;
    grid = ft_prepare_sourcemodel(cfg);
end

%----
% Create lead field = forward solution

% Only select MEG label --> TO ADAPT for other data than 4D 
chanlab = ft_channelselection('A*', Sgrad.label);
cfg = [];
cfg.grad = Sgrad;
cfg.headmodel = shell;
cfg.grid = grid;
cfg.reducerank = 2;
cfg.channel = chanlab;  
cfg.normalize = 'yes';
ldf_grid = ft_prepare_leadfield(cfg);

% Keep source label
ldf_grid.atlaslabel= subso.label;

model.subcortical.grid = ldf_grid;
model.subcortical.headmodel = shell;
    
function model_geom_fig(modl, pfig)
so = [];
so.pos = modl.grid.pos;
so.label = modl.grid.atlaslabel;

shell = modl.headmodel;

cp_fwd_model_fig(so, shell)

ht = title('Click on source to display its atlas region', 'color', [1 1 1], 'fontsize', 14);
export_fig([pfig, filesep, 'forward_model_geom.png'], '-m2');
saveas(gcf, [pfig, filesep, 'forward_model_geom.fig']) 

% Add MEG channels
delete(ht)
Sgrad = modl.grid.cfg.grad;
plot_megchan(Sgrad)
view(180, 0)

export_fig([pfig, filesep, 'forward_model_geom_megsens.png'], '-m2');
saveas(gcf, [pfig, filesep, 'forward_model_geom_megsens.fig']) 

close 