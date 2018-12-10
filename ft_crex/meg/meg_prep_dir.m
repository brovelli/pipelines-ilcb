function sdir = meg_prep_dir(opt)
% Prepare analysis directory name according to the options for preprocessing
% (filter, ica rejection and resampling)
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

if ~isempty(opt.epoched.resample_fs)
    sdir = [sdir, '_rs', num2str(opt.epoched.resample_fs)];
end