subj=4;
datainfo;
contrast = 'attended';
method = 'glm_perm'; %cosinefit_perm, glm_perm
hemis=[1];
freqs = 4:1:30;%32:2:80];
npermfiles = 10;

centerphase = [0 1/3 2/3 1 4/3 5/3]*pi-pi;
resultsdir = [projectdir 'results/collapsephasebin/'];
% resultsdir = 'P:/3011085.02/phasecode/results/collapsephasebin/';
% loading in data
for h=hemis
  if (strcmp(contrast, 'attended') || strcmp(contrast, 'unattended'))
    hemi = sprintf('hemi%d_', h);
  else
    hemi = [];
  end
  
  % load the decoding results for every frequency
  cnt=1;
  for f=freqs
    tmp{cnt} = load([resultsdir, sprintf('%s/10/sub%02d_decoding_%sf%d_%d', contrast, subj, hemi, f)]);
    tmp{cnt} = mean(mean(tmp{cnt}.accuracy,3),1);
    
    % now do the same for random phase bins
    if npermfiles==1
      tmp2{cnt,1} = load([resultsdir, sprintf('%s/10/sub%02d_decoding_rand_%sf%d_%d', contrast, subj, hemi, f)]);
      tmp2{cnt,1} = mean(mean(tmp2{cnt,1}.accuracy,4),2);
    else
    for k=1:npermfiles
      tmp2{cnt,k} = load([resultsdir, sprintf('%s/10/sub%02d_decoding_rand_%sf%d_%d%d', contrast, subj, hemi, f,k)]);
      tmp2{cnt,k} = mean(mean(tmp2{cnt,k}.accuracy,4),2);
    end
    end
    cnt=cnt+1;
  end
  if size(tmp2,2)>1
    for k=1:numel(freqs)
      dum = tmp2(k,:);
      tmp2b{k} = squeeze(cat(1,dum{:}));
    end
  else
    for k=1:numel(freqs)
      tmp2b{k} = squeeze(tmp2{k});
    end
  end
  % restructuring decoding accuracies
  acc = cat(1,tmp{:})';
  acc_rand = cat(3,tmp2b{:});
end
clear tmp*

for f=1:numel(freqs)
  % observed data
  ampl(1,f) = fit_sine_nobins(acc(:,f), centerphase);
  
  % surrogate data
  for k=1:size(acc_rand,1)
    ampl_rand(k,f) = fit_sine_nobins(squeeze(acc_rand(k,:,f))', centerphase);
  end
end

cfg=[];
cfg.tail = 1;
cfg.clustertail = 1;
cfg.clusteralpha = 0.05;
cfg.alpha = 0.05;
cfg.dim = [1 numel(freqs)];
cfg.feedback = 'yes';

cfg.numrandomization = 'all';
cfg.clusterstatistic = 'max';
cfg.clusterthreshold = 'nonparametric_common';
cfg.multivariate = 'yes';
stat = clusterstat(cfg, ampl_rand', ampl');


