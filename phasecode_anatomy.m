
datainfo;
subj=input('subject number (1-10)?');
% read the mri and save as mgz
mri = anatomy_dicom2mgz(subj);

% transform to mni space
[mri_resliced, transform_vox2mni] = anatomy_mgz2mni(subj, mri);

% transform to ctf space
[~, ~, transform_vox2ctf] = anatomy_mgz2ctf(subj, mri_resliced, transform_vox2mni);

% deface
mri_defaced = anatomy_deface(subj, mri_resliced, transform_vox2ctf);

% create headmodel
headmodel = anatomy_headmodel(subj, mri_defaced, transform_vox2ctf);


ses=2;
filename = [datadir, sprintf('sub%02d/meg%02d/sub%02d-meg%02d_cleandata.mat', subj, ses, subj, ses)];
load(filename, 'data')
figure
ft_plot_sens(data.grad, 'style', '*b');
hold on
ft_plot_headmodel(ft_convert_units(headmodel,'cm'));

% create 3d source model
sourcemodel = anatomy_sourcemodel3d(subj, mri_defaced, transform_vox2ctf);

%% 2d sourcemodel
datainfo
mri = ft_read_mri([projectdir, sprintf('results/anatomy/sub%02d/mni_resliced.mgz',subj)]);
% load([projectdir, sprintf('results/anatomy/sub%02d/transform_vox2ctf.mat', subj)]);
% mri.transform = eye(4);
% mri = ft_transform_geometry(transform_vox2ctf, mri);
anatomy_skullstrip(subj, mri);

qsubfeval(@qsub_anatomy_freesurfer,subj, 1,...
  'memreq', 1024^3 * 6, 'timreq', 720*60, 'batchid', sprintf('anatomy_freesurfer1_%d', subj))

% Check-up and white matter segmentation cleaning if needed
anatomy_volumetricQC(subj)

anatomy_wmclean(subj)

% Freesurfer qsub2
qsubfeval(@qsub_anatomy_freesurfer, subj, 2,...
  'memreq', 1024^3 * 7, 'timreq', 720*60, 'batchid', sprintf('anatomy_freesurfer2_%d', subj));

% Post-processing Freesurfer script: workbench HCP tool
!module load hcp-workbench
qsubfeval(@qsub_anatomy_freesurfer, subj, 3,...
  'memreq', 1024^3 * 6, 'timreq', 480*60, 'batchid', sprintf('anatmomy_workbench_%d', subj));

% Sourcemodel
anatomy_sourcemodel2d(subj);