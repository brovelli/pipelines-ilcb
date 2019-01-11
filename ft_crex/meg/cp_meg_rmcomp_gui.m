function Sdb = cp_meg_rmcomp_gui(Sdb)
% Mark bad sensor(s) and write the bad sensor selection in the param_txt file
% rmsens_*.txt
%
%-CREx-180726

Ns = length(Sdb);

for i = 1 : Ns
    dps = Sdb(i);
    dpmeg = dps.meg;
    Sprep = dpmeg.preproc;
    % If the initial data visualisation was already done
    if ~any(Sprep.do.rmc)
        continue;
    end
    
    idnam = dps.sinfo;
    
    drun = dpmeg.run.dir;
    Nr = length(drun);
    disp_subj(dps.info)
    
    Sdisp = [];
    %-- Interactive figure of spectra / channels to make a first selection of
    % bad channels - for each run (depending on rm_sens_run option, all bad
    % channels will be merged for all runs (case=='same')
    for j = 1 : Nr
        if ~Sprep.do.rmc(j)
            continue;
        end
        srun = drun{j};
        
        stit = [idnam,' -- ', srun];
        Spar = Sprep.param_run{j};
        
        rmc = Spar.rm.comp;
        
        pica = Spar.dir.ica;

        pfig = [pica, filesep, 'ica_fig'];     
        
        Sdisp.title = {'Bad ICA component selection' ; stit};    
    
        % Good channels
        Ncomp = Spar.Ncomp;
        Sdisp.good.clist = (1 : Ncomp)';
        
        % Bad channels
        Sdisp.bad.clist = rmc;  
        
        % Figures directory
        Sdisp.dir = pfig;
        
        uiwait(msgbox({'\fontsize{12}Please select the bad component(s) for subject: ';...
        ['\fontsize{13}\bf ', prep_tex(stit)]}, 'Bad components', 'help',...
        struct('WindowStyle', 'non-modal', 'Interpreter', 'tex')));
    
        [rmc, isbad] = preproc_select(Sdisp);              
        if isbad
            update_valrun(dpmeg.run.valtxt{j}, ~isbad);
            rmc = [];
            % Update valid run indication
            Sdb(i).meg.run.valid(j) = ~isbad;
            % Set all computations to 0
            Sprep.do = init_do(Sprep.do, j);
        end  
        
        Spar.rm.comp = rmc';
        
        % Write rmc
        write_bad(Sprep.param_txt.rmc.(srun), rmc, 'comp')

        Sprep.param_run{j} = Spar;
        Sprep.do.rmc(j) = 0;
    end
    
    % Merge the bad channel(s) already identified (== those with signal == 0)
    Sdb(i).meg.preproc = Sprep;
end
