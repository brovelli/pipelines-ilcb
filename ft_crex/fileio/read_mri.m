function mri = read_mri(pmri, mri_unit)
% Read MRIs from BrainVisa database with ft_read_mri
% + unzip gz folder for Windows platform
% + Check for the units and set them to mri_unit

% Everything in mm
if nargin < 2 || isempty(mri_unit)
    mri_unit = 'mm';
end

% Unzip volume file with gunzip Matlab function if OS is Windows
if ispc && strcmp(pmri(end-1:end), 'gz')
    ptest = gunzip(pmri);
    pmri = ptest{1};
end

% ft_tool('freesurfer', 'add');
% Home-made function use to prevent bug with the ft_hastoolbox function that seems occured
% when other toolbox hab been added in matlab path ?? (spm12 ?)
% >> Error using ft_hastoolbox (line 450) - the FREESURFER toolbox is not installed, see http://surfer.nmr.mgh.harvard.edu/fswiki
% >> Error in ft_read_mri (line 393) >> ft_hastoolbox('freesurfer', 1);

mri = ft_read_mri(pmri);

% Check for unit (cf. Fardin Afdideh recommendations)
% Convert units to mri_unit if not set yet or if units are differents
fprintf('\nSet units in %s\n', mri_unit);  
mri = ft_convert_units(mri, mri_unit);

