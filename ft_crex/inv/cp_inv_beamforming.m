function Sdb = cp_inv_beamforming(Sdb, bopt)
% Beamforming with LCMV (time domain) or DICS (time-frequency domain)
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
    
	waitbar(i/Np, wb, ['Beamforming: ', sinfo]);
    for j = 1 : Nr
        %%% TO DO Check if computation already done with the same options
        % if 
        % continue
        % end
        cond = psubj.meg.preproc.param_run{j}.conditions;
        % Check for beamforming option according to conditions to process /
        % set parameters for each condition
        opt = check_bmf_opt(bopt, cond);
        srun = rdir{j};
        waitbar(i/Np, wb, ['Leadfield: ', sinfo, '--', srun]);
        
        % Head model path (volume, grid and leadfield)
        pmod = psubj.fwd.model_run{j};
        
        % Raw MEG data directory
        praw = psubj.meg.continuous.raw{j};
        pmeg = psubj.meg.clean_mat{j};
        if isempty(pmeg)
            warning('Preprocessed MEG data required for beamforming computation...')
            warning('Computation abort for subject %s\n', [sinfo, '--', srun]);
            continue
        end
        % Preproc MEG data directory
        prep = fileparts(pmeg);


        % Load model data
        pmod = [prep, filesep, 'fwd_model.mat']; 
        save(pmod, 'fwd_model')

        % Head model (volume, grid and leadfield)
        pmod = Sdb(i).fwd.model_run{j};
        
        fwd_model = loadvar(pmod);
        
 
    end
      
end
close(wb);

% Set default sources analysis (LCMV)
function opt = check_bmf_opt(bopt, cond)
dopt = struct('meth', 'lcmv', 'lcmv', []);
% Check for conditon - should be at least some of the one in cond cell
if isempty(bopt.cond_in)
    bopt.cond_in = cond;
else
    if any(~ismember(bopt.cond_in, cond))
        error(['Condition definition for beamforming %s\n',...
            'not consistant with preprocessed MEG conditions'], 'mopt.bmf.cond_in');
    end
end
        
Nc = length(bopt.cond_in);

                     
mopt.bmf.meth = 'dics'; 
% For DICS, it is possible to do several time-frequency analysis with different
% parameters from the same epoched data (condition). In this case, it is
% mandatory to give the output condition names related to each kind of analysis (cond_out).
% Condition names of epoched data to use for source analysis 
% (if empty, the mopt.epoched.conditions or default conditions will be used)
% Condition names appear in the first level of the cleanTrials.mat epoched data
mopt.bmf.cond_in = {'S', 'A', 'R', 'SAR'};
% Condition names for the source analysis outputs
% If cond_out is empty, the same conditions will be defined as source analysis
% output.
mopt.bmf.cond_out =  [];
% Sliding windows center times vectors / conditions 
% If only one row (vector or cell), the same sliding window times will be
% used for all conditions
mopt.bmf.dics.twin_s = { -1.8  : 0.005 : 1.5    % Stim
                         -1.2  : 0.005 : 1      % Action
                        -1.3  : 0.005 : 1.5     % Reward
                        -0.2  : 0.01  : 3  };   % SAR
% It is possible to define specific parameters for each output condition
mopt.bmf.dics.freq_bd = [6 20 30 40 70 80 90 100 110 120];
mopt.bmf.dics.Df = [3 8 10 10 20 20 20 20 20 20];
mopt.bmf.dics.Dt = [0.50 0.25 0.20 0.20 0.10 0.10 0.10 0.10 0.10 0.10];

% Options for LCMV - only the cropping windows + bsl ?
mopt.bmf.lcmv.crop_win = [-0.200 1.500];
mopt.bmf.lcmv.bsl = [-0.200 0];