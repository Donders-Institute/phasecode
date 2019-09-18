function [mri] = anatomy_deface(subj, mri, transform)

datainfo
mridir = [projectdir, sprintf('results/anatomy/sub%02d/', subj)];

s = 'mri_ctf';
mri.transform = eye(4);
mri = ft_transform_geometry(transform, mri);


cfg=[];
cfg.method = 'spm';
mri = ft_defacevolume(cfg, mri);
ft_defacevolume([], mri);

save([mridir, 'mri_defaced'], 'mri');