
subj = 2;
filename = '/project/3011085.02/phasecode/DATA/sub02/sub02-mri/TOMMAR_EELSPA_SACCADE_77704.MR.TOMMAR_SKYRA.0007.0001.2014.09.11.09.38.14.448977.370847816.IMA';
mri = ft_read_mri(filename);
mri = ft_volumerealign([], mri);

cfg           = [];
cfg.output    = 'brain';
segmentedmri  = ft_volumesegment(cfg, mri);

cfg = [];
cfg.method='singleshell';
headmodel = ft_prepare_headmodel(cfg, segmentedmri);

filename = '/project/3011085.02/phasecode/DATA/sub02/sub02-meg01/sub02-meg01_cleandata.mat';
tmp = load(filename);

figure
ft_plot_sens(tmp.data.grad, 'style', '*b');

hold on
ft_plot_vol(ft_convert_units(headmodel,'cm'));