function Sdb = cp_meg_cleanup_art(Sdb, opt)
% First step for continuous MEG data preprocessing
%
%*- Automatic identification of channels with very high or low amplitude 
% (depending on the deviation from the median value in the [2 200] Hz band) 
%
%*- Semi-automatic detection of strong artefacts (which impact a large part of the sensors)
%
%-----
% A high pass filter with fc = 0.5 Hz is applied to the dataset to facilitate the visualization
% and the identification of recording issues that could be hidden by the low-frequency component.
%
%-CREx-181223

% Filtered version of data for bad channels / artefact identification
dopt = [];
dopt.type = 'hp';
dopt.fc = 0.5;   
dopt.res_fs = 400;
% Option to detect bad channel from spectra
% Factor threshold
dopt.spdetc.thr = 4;
% Frequency band to check for mean spectra value
dopt.spdetc.bd = [2 200];
dopt.rm_sens_run = 'each';

if nargin < 2
    opt = dopt;
else
    opt = check_opt(opt, dopt);
end

isa_s = strcmp(opt.rm_sens_run, 'same');

Ns = length(Sdb);

% Initialize waitbar
wb = waitbar(0, 'Strong artefact detection...', 'name', 'MEG preprocessing');
wb_custcol(wb, [0 0.6 0.8]);
for i = 1 : Ns
    dps = Sdb(i);
    dpmeg = dps.meg;
    Sprep = dpmeg.preproc;
    if ~any(Sprep.do.rma)
        continue;
    end
    
    % Subject info
    sinfo = Sdb(i).sinfo;
    
    % Run directories
    rdir = dpmeg.run.dir;
    Nr = length(rdir);

    badc = zeros(Nr, 1);
    for j = 1 : Nr
        % If the initial data visualisation was already done
        if ~Sprep.do.rma(j)
            continue;
        end
        Spar = Sprep.param_run{j};
		srun = rdir{j};
        stit = [sinfo, ' --  ', srun];
        opt.info = stit;
		
        waitbar((i-1)/Ns + (j-1)/(Nr*Ns), wb, ['Artefact detection: ', stit]);
        
        pclean = make_dir(Spar.dir.cleanup_fig);               
        
        %-- Filter (at least HP > 0.5 Hz to remove slow oscillations that can
        % make artifacts more difficult to see
        opt.figflag = 0;
        ftData = cmeg_extract_filt(Spar.dir.raw, opt);     

        if isempty(ftData)
            warning('Problem while extracting data: %s', Spar.dir.raw)
            continue;
        end 
        
        % Stacked spectrum with default paramaters (n=30 and dur=20 s)  
        spData = cmeg_fft_mean(ftData);
        
        %-- Spectrum with extremum values at frequency band [2 200] Hz and +/-
        % 4*std
        zchan = bad_spectrum(spData, ftData, opt.spdetc);
        
        %-- Detecting strong artefacts (only if no one was detected before)
        %-- Define a resampling version to reduce calculation time        
        cfg = [];
        cfg.resamplefs = opt.res_fs;
        resData = ft_resampledata(cfg, ftData); 
        
        opt.rmchan = zchan;        
        wart = cmeg_artefact_detc(resData, opt);
        write_bad(Sprep.param_txt.rma.(srun), wart, 'art');
        
        % Figures showing artefact correction 
        if ~isempty(wart)
            % Remove artefact to diplay continuous data plots (including bad
            % channels)
            resDatai = resData;
            resData = cmeg_artefact_rm(resDatai, wart);
            
            % Show some data with and without artefact padding
            pfig_art = make_dir([pclean, fsep, 'artefact']);
            opt.savepath = pfig_art;
            cmeg_artefact_fig(resData, resDatai, opt);
            
            % Keep artefact windows info
            Spar.rm.art = wart;            
        end
        
        %-- Add channel with extremum amplitude values
        if ~isempty(zchan)
            badc(j) = 1;
            rms = unique([Spar.rm.sens ; zchan]);
            Spar.rm.sens = rms;
        end
        
        Sprep.param_run{j} = Spar;
        Sprep.do.rma(j) = 0;
    end
    % Merge the bad channel(s) accross run
    if any(badc) 
        Sprep = merge_bad_sens(Sprep, isa_s, 0);
    end    
    Sdb(i).meg.preproc = Sprep;
end
close(wb);
