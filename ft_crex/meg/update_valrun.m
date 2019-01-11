function isval = update_valrun(ptxt, isval)
% Get or set the badrun statut (as isbad.mat file in the MEG info directory
if nargout
    isval = 1;
end
if isempty(ptxt)    
    return;
end
% Read the validity value
if nargin < 2 && exist(ptxt, 'file')
    isval = load(ptxt, '-ascii');
else
    % Write the validity value
    fid = fopen(ptxt, 'w');
    fprintf(fid, '%d', isval);
    fclose(fid);
end