for subj=1:10
datainfo;
contrast = 'attended';
method = 'glm_perm'; %cosinefit_perm, glm_perm
hemis=[1 2];
freqs = [4:1:30];% 32:2:80];
npermfiles = 10;

resultsdir = [projectdir 'results/collapsephasebin/'];
% resultsdir = 'P:/3011085.02/phasecode/results/collapsephasebin/';
% loading in data
hcnt=1;
for h=hemis
  if (strcmp(contrast, 'attended') || strcmp(contrast, 'unattended'))
    hemi = sprintf('hemi%d_', h);
  else
    hemi = [];
  end
  
  % load the decoding results for every frequency
  cnt=1;
  for f=freqs
    tmp{cnt} = load([resultsdir, sprintf('%s/sub%02d/sub%02d_decoding_%sf%d', contrast,subj, subj, hemi, f)]);
    tmp{cnt} = mean(mean(tmp{cnt}.accuracy,3),1);
    
    % now do the same for random phase bins
    if npermfiles==1
      tmp2{cnt,1} = load([resultsdir, sprintf('%s/sub%02d/sub%02d_decoding_rand_%sf%d_', contrast, subj, subj, hemi, f)]);
      tmp2{cnt,1} = mean(mean(tmp2{cnt,1}.accuracy,4),2);
    else
      for k=1:npermfiles
        tmp2{cnt,k} = load([resultsdir, sprintf('%s/sub%02d/sub%02d_decoding_rand_%sf%d_%d%d', contrast, subj, subj, hemi, f,k)]);
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
  acc(hcnt,:,:) = cat(1,tmp{:})';
  acc_rand(hcnt,:,:,:) = cat(3,tmp2b{:});
  hcnt=hcnt+1;
end
clear tmp*

zx = size(acc,2);
centerphase = [0:2/zx:(zx-1)/(zx/2)]*pi-pi;


cfg=[];
cfg.cosinefit.statistic = 'complex';
for h=1:hcnt-1
  % observed data
  s = statfun_cosinefit(cfg, squeeze(acc(h,:,:))', centerphase);
  amp(h,:) = abs(s.stat);
  ang(h,:) = angle(s.stat);
  
  % surrogate data
  tmp = reshape(permute(acc_rand, [1 3 4 2]), hcnt-1, length(centerphase), []);
  s = statfun_cosinefit(cfg, squeeze(tmp(h,:,:))', centerphase);
  amp_rand(h,:,:) = reshape(abs(s.stat), numel(freqs),[]);
  ang_rand(h,:,:) = reshape(angle(s.stat), numel(freqs),[]);
end

arand{subj} = amp_rand;
a{subj} = amp;
keep subj a arand a2 arand2 freqs
end
% [subj mean(mean(mean(acc)))]

% cfg=[];
% cfg.tail = 1;
% cfg.clustertail = 1;
% cfg.clusteralpha = 0.05;
% cfg.alpha = 0.05;
% cfg.dim = [1 numel(freqs)];
% cfg.feedback = 'yes';
% 
% cfg.numrandomization = 'all';
% cfg.clusterstatistic = 'maxsum';
% cfg.clusterthreshold = 'nonparametric_individual';
% cfg.multivariate = 'yes';
% for k=1:2
%   stat(k) = clusterstat(cfg, squeeze(amp_rand(k,:,:)), amp(k,:)');
% end

