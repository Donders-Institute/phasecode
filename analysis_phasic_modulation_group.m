function analysis_phasic_modulation_group(varargin)

doparc = ft_getopt(varargin, 'doparc', false);
hemis=[1 2];
datainfo

for subj=valid_subjects
  filename = [projectdir, sprintf('results/modulation/sub%02d_phasicmodulation_decoding', subj)];
  if doparc
    filename = [filename, '_parc'];
  end
  tmp = load(filename, 'amp', 'amp_rand', 'ang');
  for h=hemis
    amp{h}(:,subj) = tmp.amp(h,:)';
    amprand{h}(:,subj,:) = squeeze(tmp.amp_rand(h,:,:));
    ang{h}(:, subj) = tmp.ang(h,:)';
  end
end

nfreq = size(amp{1},1);

% do 2nd level permutation
cfg=[];
cfg.numrandomization = 1000;
cfg.uvar = 2;
cfg.ivar = 1;
cfg.tail = 1;
cfg.clustertail = 1;
cfg.clusteralpha = 0.05;
cfg.alpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.clusterthreshold = 'nonparametric_individual';
cfg.correctm = 'cluster';
cfg.connectivity = 1;
cfg.dim = [1 nfreq];

for h=hemis
  stat(h) = statfun_perm(cfg, amp{h}, amprand{h});
end

filename = [projectdir, 'results/stat_phasicmodulation_decoding'];
if doparc
  filename = [filename, '_parc'];
end
save(filename, 'stat', 'amp', 'amprand')
