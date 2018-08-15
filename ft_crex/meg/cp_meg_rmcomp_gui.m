function Sdb = cp_meg_rmcomp_gui(Sdb)
% Mark bad sensor(s) and write the bad sensor selection in the param_txt file
% rmsens_*.txt
%
%-CREx-180726

Ns = length(Sdb);

for i = 1 : Ns
    dps = Sdb(i);
    
    Sprep = dps.meg.preproc;
    % If the initial data visualisation was already done
    if ~any(Sprep.new_ica)
        continue;
    end
    
    idnam = dps.sinfo;
    
    drun = dps.meg.rundir;
    Nr = length(drun);
    disp_subj(dps.info)
    
    Sdisp = [];
    %-- Interactive figure of spectra / channels to make a first selection of
    % bad channels - for each run (depending on rm_sens_run option, all bad
    % channels will be merged for all runs (case=='same')
    for j = 1 : Nr
        if ~Sprep.new_ica(j)
            continue;
        end
        srun = drun{j};
        
        Spar = Sprep.param_run{j};
        
        rmc = Spar.rm.comp;
        
        pica = Spar.dir.ica;

        pfig = [pica, filesep, 'ica_fig'];     
        
        Sdisp.title = {'Bad ICA component selection' ; [idnam,' -- ', srun]};         
        
        % Good channels
        Ncomp = Spar.Ncomp;
        Sdisp.good.clist = (1 : Ncomp)';
        
        % Bad channels
        Sdisp.bad.clist = rmc;  
        
        % Figures directory
        Sdisp.dir = pfig;
        
        rmc = preproc_select(Sdisp);
        Spar.rm.comp = rmc';
        
        % Write rms
        write_bad(Sprep.param_txt.rmc.(srun), rmc, 'comp')
        
        Sprep.param_run{j} = Spar;
        Sprep.new_ica(j) = 0;
    end
    
    % Merge the bad channel(s) already identified (== those with signal == 0)
    Sdb(i).meg.preproc = Sprep;
end

function disp_subj(dps, rdir)
if nargin < 2
    rdir = [];
end
sep = '---------------------------------';
idnam = [dps.proj, ' ', dps.group, ' ', dps.subj, ' ', rdir];
fprintf('\n\n%s\n   %s\n%s\n\n', sep, idnam, sep);
