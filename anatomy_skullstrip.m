function anatomy_skullstrip(subj, mri)

% FSL variables
threshold       = 0.5;
T               = inv(mri.transform);
center          = round(T(1:3,4))';

% name for the temporary nifti file
datainfo;
filedir = [projectdir, sprintf('results/anatomy/sub%02d/', subj)];
temp   = fullfile(filedir, 'nifti_tmp');
filedir = [filedir, 'preproc/'];
filename = [filedir, 'skullstrip'];

%% Skullstrip via FSL

% Convert to nifti temporarily and save;
cfg = [];
cfg.filename = temp;
cfg.filetype = 'nifti';
cfg.parameter = 'anatomy';
ft_volumewrite(cfg, mri);

% Create the FSL command-string
str = ['/opt/fsl/5.0.9/bin/bet ',temp,'.nii ',temp];
str = [str,'-R -f ',num2str(threshold),' -c ', num2str(center),' -g 0 -m -v'];

% Call the FSL command-string
system(str);

% Read the FSL-based segmentation
seg  = ft_read_mri([temp,'-R.nii.gz']);
delete([temp,'.nii']);
delete([temp,'-R.nii.gz']);
delete([temp,'-R_mask.nii.gz']);

% Save the FSL-based segmentation in .mgz
cfg = [];
cfg.filename = filename;
cfg.filetype = 'mgz';
cfg.parameter = 'anatomy';
ft_volumewrite(cfg, seg);

% Check the plot already now
skullstrip = ft_read_mri([filename '.mgz']);
cfg = [];
cfg.interactive = 'yes';
ft_sourceplot(cfg, skullstrip);

end

