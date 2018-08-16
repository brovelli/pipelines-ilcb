function write_sens(ptxt, rms)
% Write sens
if isempty(rms) && ~iscell(rms)
    rms = {''};
end
Nc = length(rms);
fid = fopen(ptxt, 'w');
if ~fid
    warning('Unable to write bad sensors in %s', ptxt);
    return;
end

if Nc == 1
    fprintf(fid, '%s', rms{1});
    fclose(fid);
    return;
end
fprintf(fid, '%s\n', rms{1: end-1});
fprintf(fid, '%s', rms{end});
fclose(fid);