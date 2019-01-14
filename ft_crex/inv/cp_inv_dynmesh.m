function cp_inv_dynmesh(Spow, pfig)
%-CREx190110
if nargin < 2 || isempty(pfig)
    pfig = make_dir('dynmesh_fig', 1);
end

Spow = Spow.cortical;
% Info on frequency band / param
si = Spow.info.param;
sparam = [upper(si.name), '/', num2str(si.foi), 'Hz'];
subj = Spow.info.subj;
cond = Spow.info.cond_out;
snam = name_save(Spow.info.subj);
snam = strrep(snam, '_', '');

opt = [];  

opt.colmap = colormap_spectral;
opt.decim = 1;
opt.xlim = Spow.time([1 end]);

opt.condition = cond;
Sp = Spow;

%-- Mean across trials
Sp.pow = squeeze(mean(Spow.pow, 1)); 

opt.info = [snam, '_', 'avgtrial_', cond];
opt.title = [subj, ' - ', sparam, ' - Average across trials'];
opt.savepath = make_dir([pfig, filesep, 'avgtrials']);
cp_inv_dynmesh_fig(Sp, opt)

%-- Avg/ROI
ulab = Spow.mean_roi.label;
allab = Spow.atlas.label;
Na = length(ulab);
ins = Spow.inside;

%-- AVG across ROI and trials
[nd, nt] = size(squeeze(Spow.pow(1, :, :)));
pow = NaN(nd, nt);
for i = 1 : Na
    id = strcmp(allab, ulab{i}) & ins;
    pow(id, :) = repmat(Spow.mean_roi.mpow(i, :), sum(id), 1);
end

Sp.pow = pow;
opt.info = [snam, '_avgroi_avgtrial_', cond];
opt.title = [subj, ' - ', sparam, ' - Average across ROI and trials'];
opt.savepath = make_dir([pfig, filesep, 'avgroi_avgtrial']);
cp_inv_dynmesh_fig(Sp, opt)
