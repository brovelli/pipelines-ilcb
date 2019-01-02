function Sdb = cp_meg_cleanup_sens(Sdb, opt)
% Continuous MEG data preprocessing - selection of sensors with recording
% issues
%
%*- Figures to help for bad sensor identification
%
%-----
% A high pass filter with fc = 0.5 Hz is applied to the dataset to facilitate the visualization
% and the identification of recording issues that could be hidden by the low-frequency component.
%
% Figures are saved in cleanup_fig directory. 3 types of figures are
% generated:
%
% - fft plots - spectra are computed for all channels and superimposed on the
% same figure. This figure will be used for an interactive selection of
% spectra displaying abnormal shape.
%   [saved in 'fftplots' subdirectory]
%
% - continuous time series plots + 2 details of 10 s duration (one around the
% maximum amplitude). The data set is first resampled at 400 Hz to save time
% when printing figures
%   ['datadisp' subfolder]
%
% - the layout of MEG sensor with a color code at each channel that is an
% indicator of the signal amplitude compared to the channel that displays the
% minimum value (except for the channels with zeros values)
%   ['datadisp_filt' subfolder - with the same color code being used for the time
%   series plots]
%
%-CREx-180726

% Filtered version of data for bad channels / artefact identification
dopt = [];
dopt.type = 'hp';
dopt.fc = 0.5;   
dopt.res_fs = 400;
dopt.rm_sens_run = 'each';

if nargin < 2
    opt = dopt;
else
    opt = check_opt(opt, dopt);
end

isa_s = strcmp(opt.rm_sens_run, 'same');

Ns = length(Sdb);

% Initialize waitbar
wb = waitbar(0, 'Extract continuous data...', 'name', 'MEG preprocessing');
wb_custcol(wb, [0 0.6 0.8]);
for i = 1 : Ns
    dps = Sdb(i);
    dpmeg = dps.meg;
    Sprep = dpmeg.preproc;
    if ~any(Sprep.do.rms)
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
        if ~Sprep.do.rms(j)
            continue;
        end
        Spar = Sprep.param_run{j};
		srun = rdir{j};
        stit = [sinfo, ' --  ', srun];
        opt.info = stit;
        
        waitbar((i-1)/Ns + (j-1)/(Nr*Ns), wb, ['Cleanup figures: ', stit]);
         
        pclean = make_dir(Spar.dir.cleanup_fig);               
        
        %-- Filter (at least HP > 0.5 Hz to remove slow oscillations that can
        % make artifacts more difficult to see
        opt.figflag = 0;
        ftData = cmeg_extract_filt(Spar.dir.raw, opt);     

        if isempty(ftData)
            warning('Problem while extracting data: %s', Spar.dir.raw)
            continue;
        end  
        
        %-- Define a resampling version for continuous plots    
        cfg = [];
        cfg.resamplefs = opt.res_fs;
        resData = ft_resampledata(cfg, ftData);   
        
        % Remove the strong artefacts before making figures of continuous data
        % in time or frequency domain
        if ~isempty(Spar.rm.art)
            wart = Spar.rm.art;

            % Remove artefact to diplay continuous data plots (including bad
            % channels)
            resData = cmeg_artefact_rm(resData, wart);
         
            % Do spectra excluding bad channels and with artefact correction
            ftData = cmeg_artefact_rm(ftData, wart);                        
        end
        
        %-- Bad channels from sprectrum and data amplitude values
        spData = cmeg_fft_mean(ftData);
        bchan = unique([bad_spectrum(spData, ftData); Spar.rm.sens]);
        
        %-- Continuous data plot => based on a resampling version of the data to
        % make figure faster - all channels are kept to see the bad one too
        % before confirming the bad channel selection
        pfig_vis = make_dir([pclean, filesep, 'datadisp']); 
        Scol = cmeg_chancheck_fig(resData, pfig_vis, opt.info, bchan);
        opt.savepath = pfig_vis;
        opt.colors = Scol;
        opt.rmchan = bchan;
        cmeg_plot_data(resData, opt)

        %-- FFT plots
        opt.savepath = make_dir([pclean, filesep,'fftplots']); 
        % Figure of superimposed spectra to help for bad channel identification
        Spar.rms_fig = cmeg_fft_fig(spData, opt);
        
        fprintf('\n------\nFigures to help for channel selection:\n%s\n-------\n\n',...
            pclean);    
        
        %-- Add channel with extremum amplitude values
        if ~isempty(bchan)
            badc(j) = 1;
            Spar.rm.sens = bchan;   
            % Write rms
            write_bad(Sprep.param_txt.rms.(srun), bchan, 'sens')
        end      
        Sprep.param_run{j} = Spar;
    end
    
    % Merge the bad channel(s) accross run
    if any(badc) 
        Sprep = merge_bad_sens(Sprep, isa_s);
    end
    Sdb(i).meg.preproc = Sprep;
end
close(wb);

% Bad channel manual selection
Sdb = cp_meg_rmsens_gui(Sdb, opt.rm_sens_run);
