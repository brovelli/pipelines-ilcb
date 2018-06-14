function Sdb = cp_marsatlas(Sdb)
% Prepare atlas files (surf, vol) 
% - Realign surf and vol according to Mtrans_ref transformation matrix in the
% subject MRI_bv space
% - Realign according to fiducial previously indicated (cp_fwd_singleshell)
%
%   --> marsatlas.mat saved inside db_ft/PROJ/SUBJ/atlas
%
%-CREx180530

Np = length(Sdb);

for i = 1 : Np
    
    %-- Data paths
    
    % Surface files from BV/FS pipeline (Left and Right)
    psurf = Sdb(i).surf;
    Nh = length(psurf);
    
    % Associated texture files (Left and Right)
    ptex = Sdb(i).tex;
    
    % Volume MarsAtlas
    pvol = Sdb(i).vol{1};
    
    % Mtrans_ref to realign surf and vol  
    ptrans = Sdb(i).trans;
    
    % Mreal of MRI subject according to fiducial
    pMreal = Sdb(i).Mreal;
    
    % Conduction volume for figure of coregistration
    pshell = Sdb(i).shell;
    
    %-- Load data 
    
    % Read MarsAtlas surface files + add textures information
    masurf = cell(Nh, 1);
    for k = 1 : Nh
        psurfh = psurf{k};
        [~, snam] = fileparts(psurfh);
        % Define if Left or Right according to surf name ('Lwhite' or 'Rwhite')
        % Normally, Lwhite is in k==1 and Rwhite in k==2 as path list is
        % returned in alphabetic order - but in case of...
        if ~isempty(strfind(snam, 'Lwhite')) %#ok
            % snam = 'Lwhite';
            namh = 'surf_L';
        else
            % snam = 'Rwhite';
            namh = 'surf_R';
        end
        masurfh = read_gii(psurfh);
        % ptex list of file in the same order as psurf a priori (if not, see to
        % find the corresponding ptex depending on snam)
        masurfh.tex = read_gii(ptex{k}); 
        masurfh.name = namh;
        masurf{k} = masurfh;        
    end
    
    % Read MarsAtlas volume
    mavol = read_mars_vol(pvol);
    
    % Transformation matrices
    % to realign surf and vol in subject MRI space
    Meach = loadvar(ptrans);
    Msu2bv = Meach{1}*Meach{2};
    Mvo2bv = inv(Meach{4}*Meach{3});

    % to realign according to fiducial
    Mreal = loadvar(pMreal);
    

    % Realign MarsAtlas
    
    %-- Surface
    % Total transform
    Mt_su = Mreal*Msu2bv;
    for k = 1 : Nh
        masurf{k} = ft_transform_geometry(Mt_su, masurf{k});
    end
    
    %-- Volume
    % Total transform
    Mt_vo = Mreal*Mvo2bv;  %#ok
    % The transform field of mavol atlas will be change according to previous
    % transfom + Mt_vol
    mavol = ft_transform_geometry(Mt_vo, mavol); 
    
    % Save atlas file
    marsatlas = [];
    marsatlas.surf = masurf;
    marsatlas.vol = mavol;
    
    patlas = [make_dir([Sdb(i).dir, filesep, 'atlas']), filesep, 'marsatlas.mat'];
    save(patlas, 'marsatlas');
    
    Sdb(i).atlas = patlas;
    
    % Load volume of condution for figure
    shell = loadvar(pshell);
    
    % Figure will be saved in coreg directory
    pcor = make_dir([Sdb(i).dir, filesep, 'coreg']);
    
    % Neeed for raw MEG data file to superimpose sensors
    pmeg = Sdb(i).meg;

    coreg_fig(marsatlas, shell, pmeg, pcor);
        
end


function coreg_fig(marsatlas, shell, pmeg, pfig)

masurf = marsatlas.surf;
mavol = marsatlas.vol;

% Determine right hemisphere for one of the specific plot
ir = cellfun(@(x) strcmp(x.name, 'surf_R'), masurf);

%-- Visualize

% Head model + surfaces
plot_meshes([{shell}; masurf], {'shellcond'; masurf{1}.name;  masurf{2}.name})
view(30, 10)
title('Single-shell conduction volume and MarsAtlas surface')
axis tight
saveas(gcf, [pfig, filesep, 'marsatlas_surf_shellcond.fig'])
export_fig([pfig, filesep, 'marsatlas_surf_shellcond.png'], '-m2')
close

% Headmodel, Right MA surface + Left MA volume
plot_vol(mavol, 1:41)
plot_meshes({masurf{ir}, shell}, {masurf{ir}.name, 'shellcond'}, [0 0.9 0.7; 1 0.4 0.4], 0)
view(122, 10)
rotate3d off
title('Click on a MarsAtlas parcel to get its name', 'color', [1 1 1])
saveas(gcf, [pfig, filesep, 'marsatlas_surfR_vol_shellcond.fig'])
export_fig([pfig, filesep, 'marsatlas_surfR_vol_shellcond.png'], '-m2')
close

% Add the MEG sensors
%- Looking for the raw MEG data in pmeg directory
draw = filepath_raw(pmeg);
if ~isempty(draw)
    % Check reg MEG + vol model only
    Sgrad = ft_read_sens(draw);

    % Everything on the same figure (volumetric MarsAtlas (Right) + the Right
    % MarsAtlas surface + the singleshell head model + the MEG sensors
    plot_vol(mavol, 1:41)
    plot_meshes({masurf{ir}, shell}, {masurf{ir}.name, 'shellcond'}, [0 0.9 0.7; 1 0.4 0.4], 0)
    plot_megchan(Sgrad)
    view(160, 10)    
    title('Cliquer sur une parcelle MarsAtlas pour afficher son nom', 'color', [1 1 1])
    saveas(gcf, [pfig, filesep, 'marsatlas_surfR_vol_shellcond_megsens.fig'])
    export_fig([pfig, filesep, 'marsatlas_surfR_vol_shellcond_megsens.png'], '-m2')
    close
end