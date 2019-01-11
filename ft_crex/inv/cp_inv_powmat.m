function cp_inv_powmat(Spow, pfig)
%
%-CREx-190110
if nargin < 2 || isempty(pfig)
    pfig = make_dir('powmat_fig', 1);
end

[typ, Nty] = get_names(Spow);
for i = 1 : Nty
    ctyp = typ{i};
    cp_inv_powmat_fig(Spow.(ctyp).mean_roi, pfig)
end
