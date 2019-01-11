function disp_subj(dps, rdir)
% Disp subject info on commend windows
if nargin < 2
    rdir = [];
end
sep = '---------------------------------';
idnam = [dps.proj, ' ', dps.group, ' ', dps.subj, ' ', rdir];
fprintf('\n\n%s\n   %s\n%s\n\n', sep, idnam, sep);