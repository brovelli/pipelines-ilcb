function Sdb = cp_meg_rmsens_gui(Sdb, rm_sens_run)
% Mark bad sensor(s) and write the bad sensor selection in the param_txt file
% rmsens_*.txt
%
%-CREx-180726

isa_s = strcmp(rm_sens_run, 'same');
% Get the original MEG channel list once
hdr = loadvar([Sdb(1).meg.info{1}, filesep, 'hdr_event']);
ochan = ft_channelselection('meg', hdr.label);

Ns = length(Sdb);


for i = 1 : Ns
    dps = Sdb(i);
    dpmeg = dps.meg;
    Sprep = dpmeg.preproc;
    % If the initial data visualisation was already done
    if ~any(Sprep.do.rms)
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
        if ~Sprep.do.rms(j)
            continue;
        end
        
        srun = drun{j};
        sinfo = [idnam, ' -- ', srun];
        Spar = Sprep.param_run{j};
        
        rms = Spar.rm.sens;
        
        pfig = Spar.rms_fig;
        
        % Avoid latex interpreter bug with prep_tex
        uiwait(msgbox({'\fontsize{12}Please select the bad channel(s) for subject: ';...
            ['\fontsize{13}\bf ', prep_tex(sinfo)]}, 'Bad spectra', 'help',...
            struct('WindowStyle', 'non-modal', 'Interpreter', 'tex')));
        
        rms = cp_meg_rmsens_spectra(pfig, rms);       
        
        Sdisp.title = {'Bad channel selection' ; sinfo}; 
        
        % Good channels
        gchan = setxor(ochan, rms);
        Sdisp.good.clist = gchan;
        
        % Bad channels
        Sdisp.bad.clist = rms;  
        
        % Figures directory
        Sdisp.dir = Spar.dir.cleanup_fig;
        
        [rms, isbad] = preproc_select(Sdisp);
        if isbad
            update_valrun(dpmeg.run.valtxt{j}, ~isbad);
            rms = [];
            % Update valid run indication
            Sdb(i).meg.run.valid(j) = ~isbad; 
            % Set all computations to 0 for run j
            Sprep.do = init_do(Sprep.do, j);
        end   
            
        Spar.rm.sens = rms;
        
        % Write rms
        write_bad(Sprep.param_txt.rms.(srun), rms, 'sens')
        Sprep.param_run{j} = Spar;
        Sprep.do.rms(j) = 0;
    end
    
    % Merge the bad channel(s) accross run
    Sprep = merge_bad_sens(Sprep, isa_s);

    Sdb(i).meg.preproc = Sprep;
end
