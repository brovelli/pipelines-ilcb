function [vval, vname] = loadvar(mpath, svar)
% Load variable inside mpath MAT-file that matches with svar name
% svar can hold '*' jocker character
%
% If several data variables are found in mpath MAT-file, the data to keep is asking at the 
% command prompt
%

if nargin < 2
    svar = '*';
end
vval = [];
vname = [];

if isempty(mpath)
    return;
end

Smat = load(mpath, svar);
[vnam, Nv] = get_names(Smat);
if Nv > 1
    warning('Several variables found in data file %s', mpath)
    vnam = select_var(vnam);
end
vname = vnam{1};
vval = Smat.(vname);


function vnam = select_var(vnam)

Nv = length(vnam);
fprintf('\nChoose the variable to be used by enter the associated number\n');

for i = 1 : Nv
    fprintf('     %d -> %s\n', i, vnam{i});   
end
ivec = input('       > ');

vnam = vnam(ivec);
