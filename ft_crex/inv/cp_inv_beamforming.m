function Sdb = cp_inv_beamforming(Sdb, opt)
% Beamforming with LCMV (time domain) or DICS (time-frequency domain)
%
%-CREx180530
if nargin < 2 || isempty(opt) || ~isfield(opt, 'method')
    warning('Beamforming options not set, the beamforming will not be done...')
end
% Select subjects with available MEG data
imeg = is_meg(Sdb);
if ~any(imeg)
    return
end

Sdbm = Sdb(imeg);
Ns = length(Sdbm);

% Initialize waitbar
wb = waitbar(0, 'Beamforming computation...', 'name', 'Source analysis');
wb_custcol(wb, [0 0.6 0.8]);

bopt = [];
bopt.param = opt.(opt.method);

for i = 1 : Ns
    psubj = Sdbm(i);
    
    % Subject info
    sinfo = psubj.sinfo;
    
    rdir = psubj.meg.run.dir;
    Nr = length(rdir);
    
    for j = 1 : Nr
        
        % Path to the forward model
        pfwd = psubj.fwd.model_run{j};
        if isempty(pfwd)
            warning('Forward model not set for subject %s -- Abort beamforming', bopt.info);
            continue;
        end
        
        srun = rdir{j};
        waitbar((i-1)/Ns + (j-1)/(Nr*Ns), wb, ['Beamforming: ', sinfo, '--', srun]);
        
        % Info to be added in figure title
        bopt.info = [sinfo, '--', srun];
        
        % Available conditions in dataset
        cond = psubj.meg.preproc.param_run{j}.conditions;
        
        bopt.acond = cond;    

        bopt.fwd = pfwd;
        
        % Path to the M/EEG data
        bopt.data = psubj.meg.analysis.clean_mat{j};
        
        % Keep information if previous clean data set for source analysis is new
        bopt.new_clean = psubj.meg.analysis.new_clean(j);
        
        % Directory for results
        bopt.dir = make_dir([psubj.meg.analysis.clean_dir{j}, filesep, 'sources_', opt.method]);    
        if strcmp(opt.method, 'dics')           
            beamforming_dics(bopt);
        elseif strcmp(opt.method, 'lcmv')
            beamforming_lcmv(bopt);
        end     
    end     
end
close(wb);
% Useless for now as nothing change in Sdbm
Sdb(imeg) = Sdbm;
