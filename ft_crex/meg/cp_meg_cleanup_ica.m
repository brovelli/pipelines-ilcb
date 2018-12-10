function Sdb = cp_meg_cleanup_ica(Sdb, Nc)
% ICA for artefact correction by component rejection
%
% The ICA is processed on a resampling (at 400 Hz) and HP filtering (0.5 Hz)
% version of the continuous data set that have been cleaning from bad sensors and strong
% artefacts.
%
% Unmixing topography and topolabel are saved in icaComp_res.mat file in 
% DB_FIELDTRIP/proj/(group)/subj/meg/continuous/prep/run_*/_ica
% Figures of components are saved in the ica_fig subfolder.
%
%
%-CREx-180726

% Number of components
if nargin < 2
    Nc = 'all';
end

fopt = [];
fopt.type = 'hp';
fopt.fc = 0.5;   
fopt.res_fs = 400;
% ICA components are computed once from a data version that has been HP filetring (> 0.5 Hz),
% resampling and cleaning from bad channels and strong artefact
%%%% TO DO: add the case when new artefact / bad sensor is identified => new ICA to process 
%%%% and add figure showing the diff before_ICA / after_ICA
% The purpose being only to remove EOG and/or ECG
opt_ica = fopt;


Ns = length(Sdb);

% Initialize waitbar
wb = waitbar(0, 'ICA...', 'name', 'MEG preprocessing');
wb_custcol(wb, [0 0.6 0.8]);

for i = 1 : Ns
    Sprep = Sdb(i).meg.preproc;
	waitbar(i/Ns, wb, ['ICA: ', Sdb(i).sinfo]);  	
    
    if ~any(Sprep.do.ica)
        continue;
    end
    
    % Subject info
    sinfo = Sdb(i).sinfo;
    
    % Run directories
    rdir = Sdb(i).meg.rundir;
    Nr = length(rdir); 

    for j = 1 : Nr
        % If the initial data visualisation was already done
        if ~Sprep.do.ica(j)
            continue;
        end
        Spar = Sprep.param_run{j};
        sdir = rdir{j};
		waitbar(i/Ns, wb, ['ICA: ', sinfo, '--', sdir]); 
        
        % Main directory for figures and ICA
        pica = Spar.dir.ica;      
        
        % Check for channels to remove
        rms = Spar.rm.sens;
        chansel = chan_sel(rms);      
        
        % Check for strong artefact to remove
        wina = Spar.rm.art;
                
        %-- Filter (at least HP > 0.5 Hz to remove slow oscillations that can
        % make artifacts more difficult to see
        opt_ica.info = [sinfo, filesep, sdir];

        opt_ica.figflag = 0;
        opt_ica.channel = chansel;
        ftData = cmeg_extract_filt(Spar.dir.raw, opt_ica);     

        if isempty(ftData)
            warning('Problem while extracting data: %s', Spar.dir.raw)
            continue;
        end
        
        %-- Remove strong artefacts        
        % Eliminate artefacts defined by wina 
        ftData = cmeg_padding_artefact(ftData, wina);  
        
        %-- Define a resampling version for continuous plots and ICA computation
        % to reduce calculation time        
        cfg = [];
        cfg.resamplefs = opt_ica.res_fs;
        resData = ft_resampledata(cfg, ftData);       
        
        % Parameters for ICA calculation

        pfig = [pica, filesep, 'ica_fig'];
        if exist(pfig, 'dir')
            delete([pfig, filesep, '*.png']);
        else
            pfig = make_dir(pfig);
        end
        
        cfg = [];
        cfg.method = 'runica'; 
        cfg.demean = 'no';
        cfg.numcomponent = Nc;
        comp = ft_componentanalysis(cfg, resData);
        
        opt_ica.savepath = pfig;
        % opt.mode = mod; % => cfg for manual inspection of ICA
        cmeg_ica_topo_fig(comp, opt_ica)       
        
        % Only the comp.unmixing and topolabel are needed to reject component
        % from the original (not resampled) data set
        icaComp_res = [];
        icaComp_res.unmixing = comp.unmixing;   %%% TO DO : check if the same as pinv(comp.topo) - see meg_prep of Andrea
        icaComp_res.topolabel = comp.topolabel;
        save([pica, filesep, 'icaComp_res.mat'], 'icaComp_res');
        
        % New number of components (initialized in cp_meg_prep_init
        Spar.Ncomp = length(comp.topolabel);

        % Keep the parameters used to prepare data for ICA
        preproc_ica = [];
        % Sensor removing
        preproc_ica.rms.before = rms;
        preproc_ica.rms.after = [];
        
        % Strong artefact removing
        preproc_ica.rma.before = wina;
        preproc_ica.rma.after = [];
        preproc_ica.opt = fopt;
        save([pica, filesep, 'preproc_ica.mat'], 'preproc_ica'); 
        Sdb(i).meg.preproc.param_run{j} = Spar;     
    end  
end
close(wb);

% Input of components to reject
Sdb = cp_meg_rmcomp_gui(Sdb);
