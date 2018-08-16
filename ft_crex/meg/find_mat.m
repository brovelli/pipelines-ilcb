function [pdat, ndat, datyp] = find_mat(dpath, datopt)
% Special function to find MAT data file with the specific name :
%                  [datatype_prefix,'*',processing_suffix,'*.mat']
%
% dataopt.datatyp : string of the prefix data name to find
%
% Exemple : searching of the source*.mat data matrix obtained from data
% previously band-pass filtered between 1 and 25 Hz and resample at 200 Hz
% datopt = [];
% datopt.datatyp = 'source';
% datopt.preproc.hpfc    = 1;     % (Hz)
% datopt.preproc.lpfc    = 25;    % (Hz)
% datopt.preproc.resfs   = 200;   % (Hz)
% datopt.preproc.crop    = [0 0]; % (s)


fprintf('\n-- Search for data inside %s\n', dpath);

% Default parameters
dopt = struct('datatyp', '*',...
                    'preproc',[]);
                
if nargin <2 || isempty(datopt)==1
    datopt = dopt;
else
    datopt = check_opt(datopt, dopt);
end

if ~iscell(datopt.datatyp)
    datopt.datatyp = {datopt.datatyp};
end

% Data type (prefix of the data name to find)
prefix = datopt.datatyp;
Npre = length(prefix);

% Return preprocessing string according to dataopt.preproc
strproc = preproc_str(datopt.preproc);

% We have to check for each preprocessing string units separately
if ~isempty(strproc)
    cproc = strsplitt(strproc, '_');
else
    cproc = {''};
end

% Stop as soon as a matrix has been found with the good prefix
% Part of data matrix name to find
namat = search_name(['*' , cproc{1}, '*.mat']);
dpp = dir(fullfile(dpath, namat));

pdat = [];
ndat = [];

if isempty(dpp)
    warning('Dataset file not found with the specified option in %s', dpath);
    return;
end

% Among theses matrix, we are looking for the matrix with prefix
% matching with one prefix cell string (by the order of priority)    
dfil = [];
for i = 1 : Npre
    spref = prefix{i};
    namat = search_name([spref, '*', cproc{1}, '*.mat']);
    dfp = dir(fullfile(dpath, namat));
    if ~isempty(dfp)
        dfil = dfp;
        datyp = spref;
        break
    end   
end
  
if isempty(dfil)
    warning('Dataset file not found with the specified option in %s', dpath);
    return;
end

Nm = length(dfil);
if isempty(cproc{1}) || Nm==1
    [pdat, ndat] = dirlate(dpath, namat);
    return;
end

% Processing parameters to find in data name
pproc = {'hp'; 'lp'; 'res'; 'crop'};
Npp = length(pproc);
    
% Only keep mat that contains cproc
ikeep = zeros(Nm, 1);
for k = 1 : Nm
    [~, mnam] = fileparts(dfil(k).name);
    smat = strsplitt(mnam, '_');
    smat = smat(2:end);
    isam = intersect(smat, cproc);
    if ~isempty(isam) && length(isam)==length(cproc)
        ikeep(k) = 1;
    end
end
Nk = sum(ikeep);
if Nk <= 1
    if Nk==1
        namat = dfil(ikeep).name;
    end
    [pdat, ndat] = dirlate(dpath, namat);
    return;
end

% Try to exclude the mat file with preprocessing that was not specified in
% preproc options
dfil = dfil(ikeep);
Nd = length(dfil);
ikp = zeros(Nd, 1);
for k = 1 : Nd
    for i = 1 : Npp
        sproc = pproc{i};
        % Preprocessing found in data file name
        if containss(mnam, sproc) && ~containss(strproc, sproc)
           ikp(k) = 0;                         
        end
    end
end
Nkk = sum(ikp);
if Nkk >= 1
    dfil = dfil(ikp).name;
end
vdat = [dfil(:).datenum]';
ikp = find(vdat == max(vdat));
namat = dfil(ikp).name; %#ok
[pdat, ndat] = dirlate(dpath, namat);  


function fname = search_name(fname)

while ~strcmp(fname, strrep(fname, '**', '*'))
    fname = strrep(fname, '**', '*');
end
