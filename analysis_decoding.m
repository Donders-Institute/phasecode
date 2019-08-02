if ~exist('do_allses', 'var'), do_allses = 1; end
if ~exist('do_baselinecorrect', 'var'), do_baselinecorrect=0; end
if ~exist('do_pca', 'var'), do_pca = 0; end
if ~exist('do_correcttrials', 'var'), do_correcttrials = 0; end
if ~exist('do_phasebin', 'var'), do_phasebin = 0; end
if ~exist('do_timeresolved', 'var'), do_timeresolved = 0; end
if ~exist('contrast', 'var'), contrast = 'congruent'; end
if ~exist('do_avgtrials', 'var'), do_avgtrials = false; end
if ~exist('do_prewhiten', 'var'), do_prewhiten = false; end
if ~exist('do_smooth', 'var'), do_smooth = false; end
if ~exist('do_enet', 'var'), do_enet = false; end
if ~exist('do_cichy', 'var'), do_cichy = false; end
if ~exist('randnr' ,'var'), randnr = []; end
if ~exist('tbin','var'), tbin = randnr; end
subj=4;
rng('shuffle');

load data
datainfo;
if do_allses
  for k=1:numel(subjects(subj).sessions)
    ses = subjects(subj).sessions(k);
    filename = [datadir, sprintf('3016045.07_matves_%03d_%03d', subj, ses), '/cleandata.mat'];
    filename = [projectdir, sprintf('data/sub%02d/sub%02d-meg%02d/sub%02d-meg%02d_cleandata.mat', subj, subj, ses, subj, ses)];
    tmp{k} = load(filename, 'data');
    tmp{k} = tmp{k}.data;
  end
  data = ft_appenddata([], tmp{:});
else
  filename = [projectdir, sprintf('data/sub%02d/sub%02d-meg%02d/sub%02d-meg%02d_cleandata.mat', subj, subj, ses, subj, ses)];
  load(filename, 'data');
end
fs = data.fsample;

if do_phasebin
  divide trials into phase bins
  centerphase = [0 0.5 1 1.5]*pi;%[acos(1), acos(0), acos(-1)];
else
  centerphase = 0;
end
[trl, phase, distance, time] = analysis_alphaphase(data, centerphase);

toi = -0.1:1/data.fsample:1.2;
if do_timeresolved
  tlength = numel(toi);
else
  tlength=1;
end
data_orig=data;

% if phasebinning and changing bin per time point
if do_binpertimepoint
  trltime = tlength;
else
  trltime=1;
end
for bin = 1:size(trl,2)
  %     for itrltime = trltime
  cfg=[];
  cfg.trials = trl{itrltime,bin}; % FIXME not yet implemented to have different trials for different time points
  data = ft_selectdata(cfg, data_orig);
  
  if do_baselinecorrect
    cfg=[];
    cfg.baseline = 'yes';
    cfg.baselinewindow = [-0.1 0];
    data = ft_preprocessing(cfg, data);
  end
  
  % select time window of interest.
  cfg=[];
  cfg.toilim = [toi(1) toi(end)];
  data = ft_redefinetrial(cfg, data);
  
  if do_correcttrials
    cfg=[];
    cfg.trials = (data.trialinfo(:,7)==1) & (data.trialinfo(:,1)==data.trialinfo(:,4));
    data = ft_selectdata(cfg, data);
  end
  
  if do_pca % noise reduction
    cfg=[];
    cfg.method = 'pca';
    comp = ft_componentanalysis(cfg, data);
    tmpcomp = permute(cat(3, comp.trial{:}), [3 1 2]);
    % find which components explain 98% of the variance
    V=0;
    idx=1;
    while V<0.98
      V = V+comp.var(idx);
      idx=idx+1;
    end
    idx
    
    cfg=[];
    cfg.channel = comp.label(1:idx);
    data = ft_selectdata(cfg, comp);
  end
  if do_smooth
    dat1.trial = smoothdata(dat1.trial, 3, 'movmean', 4);
    dat2.trial = smoothdata(dat2.trial, 3, 'movmean', 4);
  end
  % split up conditions
  switch contrast
    case 'congruent'
      idx1 = 11;
      idx2 = 14;
    case 'attended'
      idx1 = [11 12; 11 13]; % CW:  first row corresponds to left, second to right hemifield
      idx2 = [13 14; 12 14];
      idx1 = idx(hemi,:);
      idx2 = idx(hemi,:);
  end
  cfg=[];
  cfg.trials = ismember(data.trialinfo(:,2), idx1);
  dat1 = ft_selectdata(cfg, data);
  cfg.trials = ismember(data.trialinfo(:,2), idx2);
  dat2 = ft_selectdata(cfg, data);
  
  %make timelock structure
  cfg=[];
  cfg.keeptrials = 'yes';
  dat1 = ft_timelockanalysis(cfg, dat1);
  dat2 = ft_timelockanalysis(cfg, dat2);
  ntrl1 = size(dat1.trial,1);
  ntrl2 = size(dat2.trial,1);
  ntrials = min([ntrl1, ntrl2]);
  

  for itrltime = trltime
    %   cfg=[];
    %   cfg.trials = trl{itrltime,bin}; % FIXME not yet implemented to have different trials for different time points
    %   data = ft_selectdata(cfg, data_orig);
  end
  if do_avgtrials
    groupsize = 5;
    ngroups = floor(ntrials/groupsize);
    
    tmpdat1 = zeros(ngroups, size(dat1.trial,2), size(dat1.trial,3));
    tmpdat2 = zeros(ngroups, size(dat2.trial,2), size(dat2.trial,3));
    idx=1;
    rand1 = randperm(ntrl1);
    rand2 = randperm(ntrl2);
    for k=1:ngroups
      tmpdat1(k,:,:) = mean(dat1.trial(rand1(idx:idx+groupsize-1),:,:),1);
      tmpdat2(k,:,:) = mean(dat2.trial(rand2(idx:idx+groupsize-1),:,:),1);
      idx=idx+groupsize;
    end
    dat1.trial = tmpdat1;
    dat2.trial = tmpdat2;
  else
    ngroups = ntrials;
  end
  
  if do_prewhiten
    cfg=[];
    cfg.covariance = 'yes';
    cfg.removemean = 'no';
    cov1 = ft_timelockanalysis(cfg, dat1);
    cov2 = ft_timelockanalysis(cfg, dat2);
    cov(:,:,1) = cov1.cov;
    cov(:,:,2) = cov2.cov;
    cov = mean(cov,3);
    cov_inv = cov^-0.5;
    
    for k=1:numel(dat1.trial,3)
      dat1.trial(:,:,k) = dat1.trial(:,:,k)*cov_inv;
      dat2.trial(:,:,k) = dat2.trial(:,:,k)*cov_inv;
    end
  end
  
  cfgs=[];
  if do_enet
    cfgs.mva = {dml.enet};
  else
    cfgs.mva=  {dml.svm};
  end
  cfgs.method = 'crossvalidate';
  cfgs.statistic = {'accuracy'};
  cfgs.type= 'nfold';
  cfgs.nfolds = size(dat1.trial,1);
  cfgs.design = [ones(size(dat1.trial,1),1); 2*ones(size(dat1.trial,1),1)];
  if do_enet
    allt=tbin;
  elseif do_binpertimepoint
    allt=1;
  else
    allt=1:tlength;
  end
  for t=allt
    if do_timeresolved
      twindow = 1/fs;
      nsample = twindow*fs;
      cfgs.latency = [toi(t) toi(t+nsample-1)];
    end
    if do_enet
      model = dml.enet('family', 'binomial', 'df', 0, 'alpha', 0.1 );
      
      for itrl=1:ngroups
        TrainData = [dat1.trial(setdiff(1:ngroups,itrl),:,t); dat2.trial(setdiff(1:ngroups,itrl),:,t)];
        TrainDesign = cfgs.design(setdiff(1:2*ngroups, [itrl itrl+ngroups]));
        model = model.train(TrainData,TrainDesign);
        TestData = [dat1.trial(itrl,:,t); dat2.trial(itrl,:,t)];
        tmp = model.test(TestData);
        accuracy(bin, itrl,t) = tmp(1,1);
        accuracy(bin, itrl+ngroups,t) = tmp(2,2);
      end
    else
      stat(bin,t) = ft_timelockstatistics(cfgs, dat1, dat2);
      accuracy(bin,t) = stat(bin,t).statistic.accuracy;
    end
    weights(bin,t) = stat(bin,t).weights;
    lambda(bin,t) = stat(bin,t).lambda;
  end
end

filename = '/project/3011085.02/phasecode/results/';
if do_phasebin
  filename = [filename, 'phasebin_svm/'];
elseif do_enet
  filename = [filename, 'enet/'];
elseif do_cichy
  filename = [filename, 'cichy/'];
else
  filename = [filename, 'svm/'];
end
filename = [filename, sprintf('sub%02d_decoding', subj)];
if exist('randnr', 'var')
  filename = [filename, sprintf('_%d', str2num(randnr))];
end

settings = struct('allses', do_allses, 'baselinecorrect', do_baselinecorrect, 'pca', do_pca, 'correcttrials', do_correcttrials, 'phasebin', do_phasebin, 'contrast', contrast, 'avgtrials', do_avgtrials, 'prewhiten', do_prewhiten);
save(filename, 'accuracy','settings')
if do_phasebin
  save(filename, 'phase', 'trl', '-append');
end