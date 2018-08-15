function Sdb = cp_meg_rmsens_gui(Sdb, copt)
% Mark bad sensor(s) and write the bad sensor selection in the param_txt file
% rmsens_*.txt
%
%-CREx-180726

isa_s = strcmp(copt.rm_sens_run, 'same');
% Get the original MEG channel list once
hdr = loadvar([Sdb(1).meg.continuous.raw{1}, filesep, 'hdr_event']);
ochan = ft_channelselection('meg', hdr.label);

Ns = length(Sdb);

for i = 1 : Ns
    dps = Sdb(i);
    
    Sprep = dps.meg.preproc;
    % If the initial data visualisation was already done
    if ~any(Sprep.new_vis)
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
        if ~Sprep.new_vis(j)
            continue;
        end
        srun = drun{j};
        Spar = Sprep.param_run{j};
        
        rms = Spar.rm.sens;
        
        pfig = Spar.rms_fig;
        rms = cp_meg_rmsens_spectra(pfig, rms);       
        
        Sdisp.title = {'Bad channel selection' ; [idnam,' -- ', srun]}; 
        
        % Good channels
        gchan = setxor(ochan, rms);
        Sdisp.good.clist = gchan;
        
        % Bad channels
        Sdisp.bad.clist = rms;  
        
        % Figures directory
        Sdisp.dir = Spar.dir.cleanup;
        
        rms = preproc_select(Sdisp);
        Spar.rm.sens = rms;
        
        % Write rms
        write_bad(Sprep.param_txt.rms.(srun), rms, 'sens')
        Sprep.param_run{j} = Spar;
        Sprep.new_vis(j) = 0;
    end
    
    % Merge the bad channel(s) already identified (== those with signal == 0)
    if isa_s
        fprintf('\nMerge the bad channels selection of the %d run(s)\n', Nr)
        rms = Sprep.param_run{1}.rm.sens;
        for j = 2 : Nr
            rms = unique([rms; Sprep.param_run{j}.rm.sens]);
        end
        write_bad(Sprep.param_txt.rms.allrun, rms, 'sens');
        
        for j = 1 : Nr
            Sprep.param_run{j}.rm.sens = rms;
        end
    end
    Sdb(i).meg.preproc = Sprep;
end


function disp_subj(dps, rdir)
if nargin < 2
    rdir = [];
end
sep = '---------------------------------';
idnam = [dps.proj, ' ', dps.group, ' ', dps.subj, ' ', rdir];
fprintf('\n\n%s\n   %s\n%s\n\n', sep, idnam, sep);
