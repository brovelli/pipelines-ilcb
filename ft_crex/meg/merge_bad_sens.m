function Sprep = merge_bad_sens(Sprep, isa_s, wflag)
if nargin < 3
    wflag = 1;
end
% Write the rmsens_allrun.txt with all bad channels toward runs
Nr = length(Sprep.param_run);
rms = Sprep.param_run{1}.rm.sens;
for j = 2 : Nr
    rms = unique([rms; Sprep.param_run{j}.rm.sens]);
end
if wflag
    write_bad(Sprep.param_txt.rms.allrun, rms, 'sens');
end

% Change the bad channel selection for each run to be the allrun selection
if isa_s
    for j = 1 : Nr
        Sprep.param_run{j}.rm.sens = rms;
    end
end