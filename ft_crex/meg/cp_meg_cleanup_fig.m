function Sdb = cp_meg_cleanup_fig(Sdb)
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

Ns = length(Sdb);

% Initialize waitbar
wb = waitbar(0, 'Extract continuous data...', 'name', 'MEG preprocessing');
wb_custcol(wb, [0 0.6 0.8]);
for i = 1 : Ns
    Sprep = Sdb(i).meg.preproc;
    if ~any(Sprep.new_vis)
        continue;
    end
    
    % Subject info
    sinfo = Sdb(i).sinfo;
	waitbar(i/Ns, wb, ['Extract: ', sinfo]);
    
    % Run directories
    rdir = Sdb(i).meg.rundir;
    Nr = length(rdir);

    for j = 1 : Nr
        % If the initial data visualisation was already done
        if ~Sprep.new_vis(j)
            continue;
        end
        Spar = Sprep.param_run{j};
		srun = rdir{j};
        opt.info = [sinfo, ' --  ', srun];
		
        waitbar(i/Ns, wb, ['Extract: ', sinfo, '--', srun]);
		
        pclean = make_dir(Spar.dir.cleanup);               
        
        %-- Filter (at least HP > 0.5 Hz to remove slow oscillations that can
        % make artifacts more difficult to see
        opt.figflag = 0;
        ftData = cmeg_extract_filt(Spar.dir.raw, opt);     

        if isempty(ftData)
            warning('Problem while extracting data: %s', Spar.dir.raw)
            continue;
        end        
        
        
        %-- FFT plots
        if ~exist(Spar.rms_fig, 'file')
            opt.savepath = make_dir([pclean, filesep,'fftplots']); 
            Spar.rms_fig = plot_spectra(ftData, opt);
        end
        
        if ~exist([pclean, filesep, 'datadisp'], 'dir')
            %-- Define a resampling version for continuous plots and ICA computation
            % to reduce calculation time        
            cfg = [];
            cfg.resamplefs = 400;
            resData = ft_resampledata(cfg, ftData);   

            %-- Continuous data plot => based on a resampling version of the data to
            % make figure faster   
            pfig_vis = make_dir([pclean, filesep, 'datadisp']); 
            [Scol, zchan] = cmeg_chancheck_fig(resData, pfig_vis, opt.info);

            %-- Add channel with zeros in bad channels selection
            Spar.rm.sens = unique([Spar.rm.sens ; zchan]);

            opt.savepath = pfig_vis;
            opt.colors = Scol;
            cmeg_plot_data(resData, opt)
        end
        
        fprintf('\n------\nFigures to help for channel selection:\n%s\n-------\n',...
            pclean);       
        
        Sdb(i).meg.preproc.param_run{j} = Spar;
    end
end
close(wb);


% Figure of superimposed spectra to help for bad channel identification
function pfig = plot_spectra(ftData, opt)

% Launch fft calculations
% Stacked spectrum with default paramaters (spsparam.n=30 and
% spsparam.dur=20 s) :
spData = meg_fftstack_calc(ftData);    

% Figures of superimposed stacked spectrum per channels
if isempty(spData.spectra)
    pfig = [];
    return;
end

% Figures with the dedicated subfunction
opt.sps_param = spData.param;
pfig = cmeg_fftstack_fig(spData, opt);
