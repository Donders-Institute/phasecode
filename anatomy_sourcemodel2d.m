function [sourcemodel] = anatomy_sourcemodel2d(subj)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

datainfo;
mridir = [projectdir, sprintf('results/anatomy/sub%02d/', subj)];
filename_sourcemodel= fullfile(mridir, 'sourcemodel2d.mat');
mridir = [mridir, 'preproc/'];
workbench_dir = fullfile(mridir, '/workbench/');


% load in the cortical sheet
filename = fullfile(workbench_dir, ['preproc', '.L.midthickness.8k_fs_LR.surf.gii']);
filename2 = strrep(filename, '.L.', '.R.');

sourcemodel = ft_read_headshape({filename, filename2});

% get the necessary coregistration information
datapath = fullfile(mridir);
load(fullfile(datapath, 'transform_vox2mni'));
T1 = transform_vox2mni;
load(fullfile(datapath, 'transform_vox2ctf'));
T2 = transform_vox2ctf;

sourcemodel = ft_transform_geometry((T2/T1), sourcemodel);
sourcemodel.inside = sourcemodel.atlasroi>0;
sourcemodel = rmfield(sourcemodel, 'atlasroi');
sourcemodel = ft_convert_units(sourcemodel, 'mm');

save(filename_sourcemodel, 'sourcemodel');

end

