function [mri, mri_ctf, transform_vox2ctf] = anatomy_mgz2ctf(subj, mri, transform)
%streams_anatomy_mgz2ctf reads in the resliced volume created with
%streams_anatomy_mgz2mni and the transformation matrix and calls
%ft_volumerealign to perform interactive coregistration to CTF coordinate
%headspace. It then saves the transformation matrix


datainfo;
mridir = [projectdir, sprintf('results/anatomy/sub%02d/', subj)];

if isequal(mri.transform, transform)
  % do nothing
else
  fprintf('adding coregistration information to the mrifile of subject %s\n', subj);
  
  mri.transform = transform;
end

%% Interactive realignment to the CTF conventions
% do an initial coregistration
cfg             = [];
cfg.interactive = 'yes';
cfg.coordsys    = 'ctf';
mri_ctf = ft_volumerealign(cfg, mri);

% save the transformation matrix
transform_vox2ctf = mri_ctf.transform;
filename_vox2ctf  = fullfile(mridir, 'transform_vox2ctf');
save(filename_vox2ctf, 'transform_vox2ctf');
  
end

