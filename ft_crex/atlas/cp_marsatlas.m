function Sdb = cp_marsatlas(Sdb)
% Prepare atlas files (surf, vol) 
% - Realign surf and vol according to Mtrans_ref transformation matrix in the
% subject MRI_bv space
% - Realign according to fiducial previously indicated (cp_fwd_singleshell)
%
%   --> marsatlas.mat saved inside db_ft/PROJ/SUBJ/atlas
%       == marsatlas structure with fields "surf" and "vol"
%
%-CREx180530

% Number of data to process
Np = length(Sdb);
% Initialize waitbar
wb = waitbar(0, 'MarsAtlas preparation...', 'name', 'Forward model');
wb_custcol(wb, [0 0.6 0.8]);
for i = 1 : Np
    waitbar((i-1)/Np, wb, ['MarsAtlas: ', Sdb(i).sinfo]);
    % Check if marsatlas already defined
    if isempty(Sdb(i).anat.atlas)
        Sdb(i) = prepare_atlas(Sdb(i));
    end     	
end
close(wb);


function dps = prepare_atlas(dps)
dpa = dps.anat;
%-- Read MarsAtlas surface files from BV/FS pipeline (Left and Right) 
% + textures information
masurf = read_mars_surf(dpa.surf, dpa.tex);
% Number of surf files
Nh = length(masurf);

%-- Read MarsAtlas volume
pvol = dpa.vol{1};
mavol = read_mars_vol(pvol);

%-- Transformation matrices
%* to realign surf and vol in subject MRI space
Meach = loadvar(dpa.trans);
Msu2bv = Meach{1}*Meach{2};
Mvo2bv = inv(Meach{4}*Meach{3});

%* to realign according to fiducial
Mreal = loadvar(dpa.Mreal);    

%-- Realign MarsAtlas

%- Surface
% Total transform
Mt_su = Mreal*Msu2bv;
for k = 1 : Nh
    masurf{k} = ft_transform_geometry(Mt_su, masurf{k});
end

%- Volume
% Total transform
Mt_vo = Mreal*Mvo2bv;  %#ok
% The transform field of mavol atlas will be change according to previous
% transfom + Mt_vol
mavol = ft_transform_geometry(Mt_vo, mavol); 
 
% Save atlas file
marsatlas = [];
marsatlas.surf = masurf;
marsatlas.vol = mavol;

patlas = [make_dir([dpa.dir, filesep, 'atlas']), filesep, 'marsatlas.mat'];
save(patlas, 'marsatlas');

dps.anat.atlas = patlas;


