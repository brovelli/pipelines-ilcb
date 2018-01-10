function bv2ft


subjects_dir = '/hpc/comco/basanisi.r/Databases/db_mne/meg_te/';
subject = 'subject_01';


% -------------------------------------------------------------------------------------------------------------------
% Anatomical files from BrainVISA
% -------------------------------------------------------------------------------------------------------------------
% White matter meshes
fname_surf_L = [ subjects_dir subject '/surf/' subject '_Lwhite.gii' ];
fname_surf_R = [ subjects_dir subject '/surf/' subject '_Rwhite.gii' ];
% MarsAtlas surface parcellation from Brainvisa
fname_tex_L = [ subjects_dir subject '/tex/' subject '_Lwhite_parcels_marsAtlas.gii' ];
fname_tex_R = [ subjects_dir subject '/tex/' subject '_Rwhite_parcels_marsAtlas.gii' ];
% MarsAtas cortical parcellatino files
fname_atlas = [ subjects_dir 'label/MarsAtlas_BV_2015.xls' ];  % Labelling xls file
fname_color = [ subjects_dir 'label/MarsAtlas.ima' ];  % Color palette
% MarsAtlas volume parcellation files
fname_vol = [ subjects_dir subject '/vol/' subject '_parcellation.nii.gz' ];
name_lobe_vol = [ 'Subcortical' ];
% Referential file list (standard files to be added)
fname_trans_ref = [ subjects_dir 'referential/referential.txt' ];
fname_trans_out = [ subjects_dir subject '/ref/' subject '-trans.trm' ];

% Read GIfTI surface file of white matter
s_L = gifti(fname_surf_L);
% Display white matter mesh
figure(1)
plot(s_L)
view(-90,180)
axis tight
lightangle(-90, 0)

% Read GIfTI surface file of MarsAtlas parcellation
t_L = gifti(fname_tex_L);
% Change label=255
t_L.cdata(t_L.cdata==255) = 0;
% Display
figure(2)
colormap lines
plot(s_L, t_L)
axis tight
view(-90,180);
pause
view(90,180);
