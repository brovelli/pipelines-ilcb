function Sdb = cp_meg_artefact(Sdb)
% Figures to help for bad channel identification
%
% A high pass filter with fc = 0.5 Hz is applied to the dataset to facilitate the visualization
% and the identification of recording issues that could be hidden by the low-frequency component.
%
% Figures are saved in _preproc_fig directory. 3 types of figures are
% generated:
%
% - fft plots - spectra are computed for all channels and superimposed on the
% same figure. 
%   [saved in 'fftplots' subdirectory]
%
% - continuous time series plots + 2 details of 10 s duration (one around the
% maximum amplitude). The data set is first resampled at 400 Hz to save time when recording figures
%   ['datadisp' subfolder]
%
% - the layout of MEG sensor with a color code at each channel that is an
% indicator of the signal amplitude compared to the channel that displays the
% minimum value (except for the channels with zeros values)
%   ['datadisp_filt' subfolder - with the same color code being used for the time
%   series plots]
%
%
%-CREx-180726

% Filtered version of data for bad channels / artefact identification
opt = [];
opt.type = 'hp';
opt.fc = 0.5;   
opt.res_fs = 400;

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
	waitbar(i/Ns, wb, ['Extract: ', sinfo]);
    
    % Run directories
    rdir = Sdb(i).meg.run.dir;
    Nr = length(rdir);

    for j = 1 : Nr
        % If the initial data visualisation was already done
        if ~Sprep.do.rma(j)
            continue;
        end
        Spar = Sprep.param_run{j};
		srun = rdir{j};
        opt.info = [sinfo, ' --  ', srun];
		
        waitbar(j/Nr, wb, ['Extract: ', sinfo, '--', srun]);
		
        pclean = make_dir(Spar.dir.cleanup_fig);               
        
        %-- Filter (at least HP > 0.5 Hz to remove slow oscillations that can
        % make artifacts more difficult to see
        opt.figflag = 0;
        opt.channel = chan_sel(Spar.rm.sens);
        ftData = cmeg_extract_filt(Spar.dir.raw, opt);     

        if isempty(ftData)
            warning('Problem while extracting data: %s', Spar.dir.raw)
            continue;
        end
        
        %-- Define a resampling version for continuous plots and ICA computation
        % to reduce calculation time        
        cfg = [];
        cfg.resamplefs = opt.res_fs;
        resData = ft_resampledata(cfg, ftData);   
         
        % Detecting strong artefacts 
        opt.rmchan = Spar.rm.sens;
        opt.force = 1;
        wart = cmeg_artefact_detc(resData, opt);
        
        write_bad(Sprep.param_txt.rma.(srun), wart, 'art');
        
        resDatai = resData;
        resData = cmeg_artefact_rm(resDatai, wart);
        % Show some data with and without artefact padding
        pfig_art = make_dir([pclean, fsep, 'artefact']);
        opt.savepath = pfig_art;
        cmeg_artefact_fig(resData, resDatai, opt);

        Spar.rm.art = wart;
        Sprep.do.rma(j) = 0;
        Sprep.param_run{j} = Spar;
    end
    Sdb(i).meg.preproc = Sprep;
end
close(wb);

% Unique artefact windows
% function art = unique_art(art)
% cpar = cellstr(num2str(art));
% cpar = cellfun(@(x) strjoint(strsplitt(x, ' '), ' - '), cpar, 'UniformOutput', 0);        
% upar = unique(cpar);
% spar = cellfun(@(x) textscan(x, '%f - %f'), upar, 'UniformOutput', 0);
% art = zeros(Na, 2);
% for i = 1 : Na
%     art(i, :) = [spar{i}{1} spar{i}{2}];
% end