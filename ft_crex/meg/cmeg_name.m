function newname = cmeg_name(ininame, suff)
% Define new matrice name according to the previous suffix indicating preprocessing indication 
%   + the new suffix to append

% Remove ext
[~, nami, ext] = fileparts(ininame);

if nargin==2 && ~isempty(suff)
    suff = name_save(suff);
    csuff = strsplitt(suff, '_')';
else
    newname = ininame;
    return;
end

snam = strsplitt(nami, '_');
Ns = length(snam);

% Core name (data type)
cnam = snam{1};

if Ns > 1
    pcsuff = snam(2:end)';
else
    newname = [cnam, '_', strjoin(csuff, '_'), ext];
    return;
end

% Remove previous processing if new one was done
sprep = {'hp', 'lp', 're'};

cp = char(pcsuff);
pcs = cellstr(cp(:, 1:2));

ipp = ismember(pcs, sprep);
prevp = pcs(ipp);

cc = char(csuff);
newcs = cellstr(cc(:, 1:2));

ipn = ismember(newcs, sprep);
newp = newcs(ipn);

if ~isempty(prevp) && ~isempty(newp)
    isk = ~ismember(prevp, newp);
    ipk = ipp==0;
    ipk(ipp) = isk;
    pcsuff = pcsuff(ipk);
end

if isempty(pcsuff)
    newname = [cnam, '_', strjoin(csuff, '_'), ext];
    return;
end

% Check for other clean-up processing [rt, rs, rc, ra]
sprep = {'rt', 'rs', 'rc', 'ra'};

Ns = length(sprep);
for i = 1 : Ns
    rr = sprep{i};
    
    isp = containss(pcsuff, rr);
    isn = containss(csuff, rr);
    if any(isp) && any(isn)
        pnum = str2double(strrep(pcsuff(isp), rr, ''));
        nnum = str2double(strrep(csuff(isn), rr, ''));
        newnum = [num2str(pnum + nnum), rr];
        pcsuff(isp) = {newnum};
        csuff(isn) = {newnum};
    end
end

rsuff = unique([csuff ; pcsuff]);

newname = [cnam, '_', strjoin(rsuff, '_'), ext];


