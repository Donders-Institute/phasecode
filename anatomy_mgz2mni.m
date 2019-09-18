function [mri_resliced transform_vox2mni] = anatomy_mgz2mni(subj, mri)
% streams_anatomy_mgz2mni reads in the files created by
% streams_anatomy_dicom2mgz, and creates the transformation matrix for conversion to MNI space
% perfroms reslicing to 256x256x256 space. It saves the reliced volume and the transformation matrix. 

%% Initialize the variables
datainfo;
mridir = [projectdir, sprintf('results/anatomy/sub%02d/', subj)];
resliced_filename   = fullfile(mridir, 'mni_resliced.mgz');


%% Transform to MNI and reslice to freesurfer-friendly dimensions

% realign to MNI space
cfg                 = [];
cfg.coordsys        = 'spm';
cfg.parameter       = 'anatomy';
cfg.method          = 'interactive';
mri_mni             = ft_volumerealign(cfg, mri);

% reslice & save the transformation matrix to the anatomy_dir
cfg                 = [];
cfg.resolution      = 1;
cfg.dim             = [256 256 256];
mri_resliced        = ft_volumereslice(cfg, mri_mni);

% Save the resliced mni-transformed mri image
cfg                 = [];
cfg.filename        = resliced_filename;
cfg.filetype        = 'mgz';
cfg.parameter       = 'anatomy';
ft_volumewrite(cfg, mri_resliced)

% Save the transformation matrix
transform_vox2mni   = mri_resliced.transform;
filename_vox2mni    = fullfile(mridir, 'transform_vox2mni');
save(filename_vox2mni, 'transform_vox2mni');

end

