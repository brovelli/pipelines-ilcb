function bv2ft


subjects_sir = '/hpc/comco/basanisi.r/Databases/db_mne/meg_te/';
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
fname_vol = subjects_dir subject '/vol/{0}_parcellation.nii.gz' ];
name_lobe_vol = [ 'Subcortical' ];
% Referential file list (standard files to be added)
fname_trans_ref = [ subjects_dir 'referential/referential.txt' ];
fname_trans_out = [ subjects_dir + '{0}/ref/{0}-trans.trm' ];