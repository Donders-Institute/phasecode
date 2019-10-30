subj=4;
datainfo;
contrast = 'attended';
method = 'glm_perm'; %cosinefit_perm, glm_perm
hemis=[1];
freqs = 4:1:30;%32:2:80];

centerphase = [0 1/3 2/3 1 4/3 5/3]*pi-pi;
resultsdir = [projectdir 'results/collapsephasebin/'];

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
    tmp{cnt} = load([resultsdir, sprintf('%s/sub%02d_decoding_%sf%d_%d', contrast, subj, hemi, f)]);
    tmp{cnt} = mean(tmp{cnt}.accuracy,3);
    
    % now do the same for random phase bins
    tmp2{cnt} = load([resultsdir, sprintf('%s/sub%02d_decoding_rand_%sf%d_%d', contrast, subj, hemi, f)]);
    tmp2{cnt} = mean(tmp2{cnt}.accuracy,4);
    cnt=cnt+1;
  end
  
  % restructuring decoding accuracies
  acc = squeeze(mean(cat(3,tmp{:}),1));
  acc_rand = squeeze(mean(cat(4,tmp2{:}),2));
end

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
cfg.clusterstatistic = 'max';
cfg.clusterthreshold = 'nonparametric_individual';
cfg.clusteralpha = 0.05;
cfg.dim = [];
stat = clusterstat(cfg, ampl_rand, ampl);


