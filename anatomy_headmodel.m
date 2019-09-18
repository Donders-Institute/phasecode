function headmodel = anatomy_headmodel(subj, mri, transform)

datainfo
mridir = [projectdir, sprintf('results/anatomy/sub%02d/', subj)];

mri.coordsys                = 'ctf';
if ~isequal(transform, mri.transform)
  warning('mri.transform was not equal to transform in function input and is replaced.')
  mri.transform               = transform;
end

cfg = [];
cfg.output = 'brain';
seg = ft_volumesegment(cfg, mri);

cfg = [];
cfg.method = 'projectmesh';
cfg.numvertices = 10000;
bnd = ft_prepare_mesh(cfg, seg);

cfg = [];
cfg.method = 'singleshell';
headmodel = ft_prepare_headmodel(cfg, bnd);
headmodel = ft_convert_units(headmodel, 'mm');

save([mridir 'headmodel.mat'], 'headmodel');

end

