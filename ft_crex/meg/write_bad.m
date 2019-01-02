function write_bad(ptxt, rmb, styp)
if nargin < 3 || isempty(styp)
    if iscell(rmb)
        styp = 'sens';
    elseif size(rmb, 2)==2
        styp = 'art';
    else
        styp = 'elements';
    end
end
if strcmp(styp, 'sens') || strcmp(styp, 'sensors')
    write_sens(ptxt, rmb);
elseif strcmp(styp, 'art')
    write_art(ptxt, rmb);
else
    write_comp(ptxt, rmb, styp);
end

    
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

% Write artefact windows
function write_art(ptxt, rma)
if isempty(rma)
    rma = [];
    Nc = 0;
else
    Nc = length(rma(:, 1));
end
fid = fopen(ptxt, 'w');
if ~fid
    warning('Unable to write strong artefact(s) in %s', ptxt);
    return;
end

for i = 1 : Nc
    fprintf(fid, '%4.4f\t%4.4f', rma(i, 1), rma(i, 2));  
    if i < Nc
        fprintf(fid, '\n');
    end
end
fclose(fid);

% Write comp/trials
function write_comp(ptxt, rmc, styp)
if strcmp(styp, 'comp')
    styp = 'components';
end
if isempty(rmc)
    rmc = [];
end
Nc = length(rmc);
fid = fopen(ptxt, 'w');
if ~fid
    warning('Unable to write bad %s in %s', styp, ptxt);
    return;
end

if Nc == 1
    fprintf(fid, '%d', rmc);
    fclose(fid);
    return;
elseif Nc > 1
    fprintf(fid, '%d\n', rmc(1:end-1));
    fprintf(fid, '%d', rmc(end));
end
fclose(fid);