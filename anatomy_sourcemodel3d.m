function [sourcemodel] = anatomy_sourcemodel3d(subj, mri, transform)
% computes a 3D regular grid with 
% specified resolution based on an inverse warp of a template
% grid in MNI space.
%

datainfo
mridir = [projectdir, sprintf('results/anatomy/sub%02d/', subj)];

resolution = 6;
fname   = fullfile('/project/3011085.02/scripts/fieldtrip/', 'template', 'sourcemodel', ['standard_sourcemodel3d',num2str(resolution),'mm.mat']);
load(fname);

mri.coordsys                = 'ctf';
if ~isequal(transform, mri.transform)
  warning('mri.transform was not equal to transform in function input and is replaced.')
  mri.transform               = transform;
end

% create the grid
cfg = [];
cfg.grid.warpmni    = 'yes';
cfg.grid.template   = sourcemodel;
cfg.grid.nonlinear  = 'yes';
cfg.mri = mri; 
sourcemodel = ft_prepare_sourcemodel(cfg);
sourcemodel = ft_convert_units(sourcemodel, 'mm');

% remove the mri-structure from grid.cfg
sourcemodel.cfg = rmfield(sourcemodel.cfg, 'mri');
sourcemodel.cfg = rmfield(sourcemodel.cfg, 'callinfo');

save([mridir, 'sourcemodel3d'], 'sourcemodel');
end