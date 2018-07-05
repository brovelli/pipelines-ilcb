function ptex = cp_hiphop_tex(pmesh)
% Determine texture file paths associated with the decimated hiphop meshes (Left and Right) 
%
% According to the number of nodes of the mesh file, the associated hiphop
% texture files paths inside the ft_crex/atlas/hiphop138/decim/*K directory are
% returned
%
% --- Input
% pmesh = path of one of the subject's mesh (left or white) - expected to be a
% GII file from brainvisa processing (subject*_*white_remeshed_hiphop.gii)
%
% --- Output
% ptex: paths to the left and right texture files stored inside ft_crex toolbox  ([2 x 1] cell)
%
%
%-CREx-180629
%

% Read the mesh file
Smesh = read_gii(pmesh);
Npt = length(Smesh.pnt(:, 1));

ndec = [1000 2000 3000 4000];
sdec = {'1K', '2K', '3K', '4K'};

dec = sdec{find(ndec <= Npt, 1, 'last')};
% or index of min distance: [~, i] = min(abs(ndec-Npt))

% Be sure to find the decimated texture file in ft_crex toolbox
ptex = fullfile(ptool, 'atlas', 'hiphop138', 'decim', dec, ['hiphop138_*white_dec_', dec, '_parcels_marsAtlas.gii']);