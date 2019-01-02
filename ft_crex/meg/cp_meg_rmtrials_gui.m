function Sdb = cp_meg_rmtrials_gui(Sdb)
% Select bad trials and write the selection in the param_txt file
% rmtrial_*.txt
%
%-CREx-180726

Ns = length(Sdb);

for i = 1 : Ns
    dps = Sdb(i);
    
    Sprep = dps.meg.preproc;
    % If the initial data visualisation was already done
    if ~any(Sprep.do.rmt)
        continue;
    end
    
    idnam = fullfile(dps.proj, dps.group, dps.subj);
    
    drun = dps.meg.run.dir;
    Nr = length(drun);
    
    Sdisp = [];
    %-- Interactive figure of spectra / channels to make a first selection of
    % bad channels - for each run (depending on rm_sens_run option, all bad
    % channels will be merged for all runs (case=='same')
    for j = 1 : Nr
        if ~Sprep.do.rmt(j)
            continue;
        end
        srun = drun{j};
        
        Spar = Sprep.param_run{j};
        
        Str = Spar.rm.trials;
        
        [cond, Nc] = get_names(Str);
        
        uiwait(msgbox({'\fontsize{12}Please select the bad trials(s) for subject: ';...
        ['\fontsize{13}\bf ', prep_tex([idnam, ' -- ', srun])]}, 'Bad trials', 'help',...
        struct('WindowStyle', 'non-modal', 'Interpreter', 'tex')));
        
        for k = 1 : Nc
            cnam = cond{k};
            Ntr = Spar.Ntr.(cnam);
            if strcmp(cnam, 'allcond')                
                cnamt = 'same for all conditions';
            else
                cnamt = ['Condition: ', cnam];
            end
            Sdisp.title = {'Bad trials selection' ; [idnam,' -- ', srun, ' / ', cnamt]};         

            % Good trials

            Sdisp.good.clist = (1 : Ntr)';

            % Bad channels
            Sdisp.bad.clist = Str.(cnam);  

            % Figures directory
            Sdisp.dir = [];
        
            rmt = preproc_select(Sdisp);
            Spar.rm.trials.(cnam) = rmt';

            % Write rms
            write_bad(Sprep.param_txt.rmt.(srun).(cnam), rmt, 'trials')
        end
        
        Sprep.param_run{j} = Spar;
        Sprep.do.rmt(j) = 0;
    end
    
    % Merge the bad channel(s) already identified (== those with signal == 0)
    Sdb(i).meg.preproc = Sprep;
end
