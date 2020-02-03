% test whether we can decode the orientation of the (un)attended from the
% MEG signal, at the group level.

datainfo;
contrasts = {'attended', 'unattended'};
filedir = [projectdir 'results/collapsephasebin/'];
hemis = [1 2];
nperm=100;


acc = rand(numel(contrasts), numel(hemis), numel(valid_subjects), nperm);
accrand = rand(numel(contrasts), numel(hemis), numel(valid_subjects), nperm);

for c=1:numel(contrasts)
  for h=hemis
    for subj=valid_subjects
      filename = fullfile([filedir, contrasts{c},'/', sprintf('sub%02d/sub%02d_decoding_hemi%d.mat', subj, subj, h)]);
      tmp = load(filename);
      acc(c,h,subj,:) = mean(tmp.accuracy,2); % mean accuracy over permutations
      
      filename = fullfile([filedir, contrasts{c},'/', sprintf('sub%02d/sub%02d_decoding_rand_hemi%d.mat', subj, subj, h)]);
      tmp = load(filename);
      accrand(c,h,subj,:) = mean(tmp.accuracy,2);
    end
  end
end


cfgs = [];
cfgs.method = 'analytic';
cfgs.statistic = 'indepsamplesT';
cfgs.design = [ones(1,numel(valid_subjects)), 2*ones(1,numel(valid_subjects))];

for c=1:numel(contrasts)
  for h=1:numel(hemis)
    dat(c,h) = [];
    dat(c,h).dimord = 'rpt_chan_time';
    dat(c,h).label = {'accuracy'};
    dat(c,h).time=0;
    dat(c,h).trial = mean(acc(c,h,:,:),4);
    
    datrand(c,h) = dat(c,h);
    datrand(c,h).trial = mean(accrand(c,h,:,:),4);
    
    stat(c,h) = ft_timelockstatistics(cfgs, dat(c,h), datrand(c,h));
  end
end



