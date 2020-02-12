% test whether we can decode the orientation of the (un)attended from the
% MEG signal, at the group level.

datainfo;
contrasts = {'attended', 'unattended'};
filedir_orig = [projectdir 'results/collapsephasebin/'];
hemis = [1 2];
nperm=100;
if ~exist('doeye','var'), doeye=false; end


acc = rand(numel(contrasts), numel(hemis), numel(valid_subjects), nperm);
accrand = rand(numel(contrasts), numel(hemis), numel(valid_subjects), nperm);

for c=1:numel(contrasts)
  for h=hemis
    for subj=valid_subjects
      filedir = fullfile([filedir_orig, contrasts{c},'/', sprintf('sub%02d/', subj)]);
      if doeye
        filedir = [filedir ,'eye/'];
      end
      filename = [filedir, sprintf('sub%02d_decoding_hemi%d.mat', subj, h)];
      tmp = load(filename);
      acc(c,h,subj,:) = mean(tmp.accuracy,2); % mean accuracy over permutations
      
      filename = [filedir, sprintf('sub%02d_decoding_rand_hemi%d.mat', subj, h)];
      tmp = load(filename);
      accrand(c,h,subj,:) = mean(tmp.accuracy,2);
    end
  end
end


cfgs = [];
cfgs.method = 'analytic';
cfgs.statistic = 'depsamplesT';
cfgs.design = [ones(1,numel(valid_subjects)), 2*ones(1,numel(valid_subjects)); 1:numel(valid_subjects), 1:numel(valid_subjects)];
cfgs.parameter = 'trial';
cfgs.alpha = 0.05/sum([numel(hemis), numel(contrasts)]);
cfgs.ivar = 1;
cfgs.uvar = 2;

for c=1:numel(contrasts)
  for h=1:numel(hemis)
    dat(c,h).dimord = 'rpt_chan_time';
    dat(c,h).label = {'accuracy'};
    dat(c,h).time=0;
    dat(c,h).trial = squeeze(mean(acc(c,h,:,:),4));
    
    datrand(c,h) = dat(c,h);
    datrand(c,h).trial = squeeze(mean(accrand(c,h,:,:),4));
    
    stat(c,h) = ft_timelockstatistics(cfgs, dat(c,h), datrand(c,h));
    accuracy(c,h).mean = mean(dat(c,h).trial);
    accuracy(c,h).std = std(dat(c,h).trial);
  end
end

filename = [projectdir, 'stat_decoding']; 
if doeye
  filename = [filename, '_eye'];
end
save(filename, 'stat', 'accuracy')



