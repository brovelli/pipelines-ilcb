function ftData = prepare_data(Spar, fopt)
% Extract continuous raw data and prepare data according to Spar.rm options
% and to fopt (filtering options)
rmc = Spar.rm.comp;
if ~isempty(rmc)
   ftData = prep_data_ica(Spar, fopt);
else
   ftData = prep_data(Spar, fopt, Spar.rm.sens, Spar.rm.art);
end
function ftData = prep_data(Spar, fopt, rms, rma)

%-- Filter
fopt.figflag = 0;
fopt.channel = chan_sel(rms);
ftData = cmeg_extract_filt(Spar.dir.raw, fopt);  

%-- Remove strong artefacts        
ftData = cmeg_artefact_rm(ftData, rma);  

function ftData = prep_data_ica(Spar, fopt)

comp = loadvar([Spar.dir.ica, filesep, 'icaComp_res.mat']);
prepi = loadvar([Spar.dir.ica, filesep, 'preproc_ica.mat']);
fopti = prepi.opt;
ftData = prep_data(Spar, fopti, prepi.rms.before, prepi.rma.before);

% Associated components
cfg            = [];
cfg.method     = 'runica';
cfg.unmixing   = comp.unmixing; %%% TO DO - check pinv(comp.topo);
cfg.topolabel  = comp.topolabel;
comp      = ft_componentanalysis(cfg, ftData); 
%%% TO DO - see if same result as comp computed in cp_meg_cleanup_ica

% ICA correction
cfg = [];
cfg.component = Spar.rm.comp;
ftData = ft_rejectcomponent(cfg, comp, ftData);

% Apply the param_after
chansel = chan_sel(prepi.rms.after);
if ~isempty(chansel)
    cfg = [];
    cfg.channel = chansel;
    ftData = ft_preprocessing(cfg, ftData);
end

%-- Remove strong artefacts        
ftData = cmeg_artefact_rm(ftData, prepi.rma.after); 

%--- Final filter if fopti ~=fopt
if ~strcmp(fopti.type, fopt.type) || ~all(fopti.fc == fopt.fc)
    ftData = cmeg_filt(ftData, fopt);
end