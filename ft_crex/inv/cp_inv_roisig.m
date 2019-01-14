function cp_inv_roisig(Spow, pfig)
%-CREx190110
if nargin < 2 || isempty(pfig)
    pfig = make_dir('roisig_fig', 1);
end

% Initialize waitbar
Swb = [];
Swb.wb = waitbar(0, 'ROI signal figures...', 'name', 'ROI-signal figures');
wb_custcol(Swb.wb, [0 0.05 0.93]);

[typ, Nty] = get_names(Spow);
for i = 1 : Nty
    ctyp = typ{i};
    Smean = Spow.(ctyp).mean_roi;
    Swb.Nt = Nty;
    Swb.k = i;
    roisig_fig(Smean, pfig, Swb);
end
close(Swb.wb);
function roisig_fig(Smean, pfig, Swb)

%- All ROI labels
allab = Smean.atlas.label;
ulab = Smean.label;
Na = length(ulab);

cola = color_group(Na);

msh = Smean.mesh;
% No face for subcortical => display the headmodel
if ~isfield(msh, 'faces')
    msh = msh.hdm;
    msh.facecolor = [0.89 1 0.89];
    msh.facealpha = 0.12;
else
    msh.facecolor = [0.3 0.85 0.8];
    msh.facealpha = 0.13;
end

myabs = [min(Smean.mpow(:)) max(Smean.mpow(:))];

time = Smean.time;

Sinfo = Smean.info;
fnam = upper(Sinfo.param.name);
cond = Sinfo.cond_out;
snam = name_save(Sinfo.subj);
snam = strrep(snam, '_', '');
subj = Sinfo.subj;

opt = [];
opt.name = [];
opt.meshes = {msh};

for i = 1 : Na 
    lab = ulab{i};
    flab = Smean.labfull{i};   
    ind = strcmp(allab, lab);      
    
    waitbar((Swb.k - 1)/Swb.Nt + (i-1)/(Swb.Nt*Na), Swb.wb, [subj, ': ', lab]);

    opt.dip_pos = Smean.mesh.vertices(ind, :);
    if isfield(msh, 'norm')
        opt.dip_ori = msh.norm(ind, :);
    end
    opt.dip_col = cola(i, :);
    opt.side = lab(end);
    
    opt.title = {[subj, ' - Region: ', flab, ' (Ndip = ',...
        num2str(Smean.Ndip(i)), ')']; ['Condition: ', cond]};    
    %-- Mean
    csig = Smean.mpow(i, :);
    if isfield(Smean, 'ci_mpow')
        opt.ci_sig = Smean.ci_mpow(i, :);  
    end
    opt.sig_title = [fnam, ' - Mean power across dipoles and trials'];
    opt.info = [snam, '_', lab, '_avgtrial_', cond];  
    opt.ylim_abs = myabs;
    
    opt.savepath = make_dir([pfig, filesep, 'avgtrials']);
    
    cp_inv_roisig_fig(time, csig, opt)
end