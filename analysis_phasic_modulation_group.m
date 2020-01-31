function analysis_phasic_modulation_group(varargin)

doparc = ft_getopt(varargin, 'doparc', 'false');

datainfo

for k=valid_subjects
  filename = [projectdir, sprintf('results/modulation/sub%02d_phasicmodulation_decoding', subj)];
  if doparc
    filename = [filename, '_parc'];
  end
  tmp{k} = load(filename, 'amp', 'amprand', 'ang');
end

%FIXME: From here on, do statistics (yet to be implemented)

cfg=[];
cfg.tail = 1;
cfg.clustertail = 1;
cfg.clusteralpha = 0.05;
cfg.alpha = 0.05;
cfg.dim = [1 numel(freqs)];
cfg.feedback = 'yes';

cfg.numrandomization = 'all';
cfg.clusterstatistic = 'maxsum';
cfg.clusterthreshold = 'nonparametric_individual';
cfg.multivariate = 'yes';
for k=1:2
  stat(k) = clusterstat(cfg, squeeze(amp_rand(k,:,:)), amp(k,:)');
end

filename = [projectdir, 'results/stat_phasicmodulation_decoding'];
if doparc
  filename = [filename, '_parc'];
end
save(filename, 'stat')
