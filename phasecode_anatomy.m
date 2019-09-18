
datainfo;
% read the mri and save as mgz
mri = anatomy_dicom2mgz(subj);

% transform to mni space
[mri_resliced, transform_vox2mni] = anatomy_mgz2mni(subj, mri);

% transform to ctf space
[~, ~, transform_vox2ctf] = anatomy_mgz2ctf(subj, mri, transform_vox2mni);

% deface
mri_defaced = anatomy_deface(subj, mri_resliced, transform_vox2ctf);

% create headmodel
headmodel = anatomy_headmodel(subj, mri_defaced, transform_vox2ctf);

% create 3d source model
sourcemodel = anatomy_sourcemodel3d(subj, mri_defaced, transform_vox2ctf);