function sdir = meg_prep_dir(opt)

fopt = opt.continuous.filt;
if ~isempty(fopt.type) && ~isempty(fopt.fc)
    fc = fopt.fc;
    if strcmp(fopt.type, 'bp')
        sfc = [num2str(fc(1)), '_', num2str(fc(2))];
    else
        sfc = num2str(fc(1));
    end
    sfc(sfc=='.') = '';
    sdir = ['filt_', fopt.type, '_', sfc];
else
    sdir = 'no_filt';
end

if opt.continuous.ica.reject
    sdir = [sdir, '_ica'];
end