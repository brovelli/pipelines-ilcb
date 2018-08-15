function Sdat = add_info(Sdat, fnam, newval)
if ~isfield(Sdat, 'hdr')
    Sdat.hdr = [];
end
if ~isfield(Sdat.hdr, 'info')
    Sdat.hdr.info = [];
end
if ~isfield(Sdat, 'preproc')
    Sdat.hdr.info.preproc = [];
end

Sprep = Sdat.hdr.info.preproc;

if ~isfield(Sprep, fnam)
    Sprep.(fnam) = newval;
else
    % Store new preproc as new cell
    vali = Sprep.(fnam);
    if ~iscell(vali)
        vali = {vali};
    end

    if ~iscell(newval)
        newval = {newval};
    end
    szi = size(vali);    
    % Add as column vector
    if szi(2)==1
        vali = vali';
    end
    if szn(1) > 1
        newval = newval';
    end
    valf = [vali, newval];        

    Sprep.(fnam) = valf;
end
Sdat.hdr.info.preproc = Sprep;
  

    