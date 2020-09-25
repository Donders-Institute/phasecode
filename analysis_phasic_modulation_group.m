function analysis_phasic_modulation_group

hemis=[1 2];
datainfo

for subj=valid_subjects
  filename = [projectdir, sprintf('results/modulation/sub%02d_phasicmodulation_decoding', subj)];
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
save(filename, 'stat', 'amp', 'amprand')

% calculate bayes factor
n=10;
cfg = [];
cfg.method = 'analytic';
cfg.statistic = 'ft_statfun_bayesfactor';
cfg.ivar = 1;
cfg.uvar = 2;
cfg.design = [ones(1,n), 2*ones(1,n); 1:n, 1:n];
for h=hemis
    dat{h} = [];
    dat{h}.label={'modulation'};
    dat{h}.dimord = 'rpt_chan_freq';
    dat{h}.powspctrm(:,1,:) = amp{h}';
    dat{h}.freq = 4:30;
    
    datrand{h} = removefields(dat{h}, 'powspctrm');
    datrand{h}.powspctrm(:,1,:) = transpose(nanmean(amprand{h},3));
    bayesfactor(h) = ft_freqstatistics(cfg, dat{h}, datrand{h});
end
save(filename, 'bayesfactor', '-append')
