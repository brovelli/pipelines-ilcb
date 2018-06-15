function pmri_filt = cp_mri_denoise(pmri)
% Launch the SPM12 process to denoise the MRI images by applying the
% 'Spatially adaptive non-local means denoising filter' (SANLM tool from CAT tool)
% Return the file path of the denoised MRI file
%
%-CREx-180411

% Add fieldtrip/external SPM toolbox 
ft_tool('spm12', 'add', 1);

% Initialise SPM
spm('defaults','fmri');  
spm_jobman('initcfg');

% Prefix to add to denoised mri file 
prfx = 'sanlm_';

% Denoising using CAT12
proc = [];
proc.spm.tools.cat.tools.sanlm.data = {pmri};
proc.spm.tools.cat.tools.sanlm.prefix = prfx;
proc.spm.tools.cat.tools.sanlm.NCstr = -Inf;
proc.spm.tools.cat.tools.sanlm.rician = 0;

[pdir, nmri, ext] = fileparts(pmri);
pmri_filt = [pdir, filesep, prfx, nmri, ext];

spm_jobman('initcfg');
spm_jobman('run', {proc});