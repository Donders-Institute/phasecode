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

% load data
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
  % divide trials into phase bins
  centerphase = [0 0.5 1 1.5]*pi;%[acos(1), acos(0), acos(-1)];
else
  centerphase = 0;
end
[trl, phase, distance, time] = analysis_alphaphase(data, centerphase);

toi = -0.1:1/data.fsample:1.2;
if do_timeresolved
  tlength = numel(toi);
  twindow = 1/fs;
  nsample = twindow*fs;
else
  tlength=1;
end

%% Data selection
% baseline correction
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

% only select validly cued and correct trials
if do_correcttrials
  cfg=[];
  cfg.trials = (data.trialinfo(:,7)==1) & (data.trialinfo(:,1)==data.trialinfo(:,4));
  data = ft_selectdata(cfg, data);
end
data_orig=data;

% If phasebinning and changing bin per time point, loop over time points.
% If phasebinning and not changing bin per time point, make all bins
% according to first time point and don't loop (loop over time will start
% later on).
if do_binpertimepoint
  trltime = 1:tlength;
else
  trltime=1;
end

%% loop over bins (only >1 when phasebinning)
for bin = 1:size(trl,2)
  % select trial subset (when phasebinning)
  for itrltime = trltime
    cfg=[];
    cfg.trials = trl{itrltime,bin}; % FIXME not yet implemented to have different trials for different time points
    data = ft_selectdata(cfg, data_orig);
    
    %% Feature reduction (PCA)
    if do_pca
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
    
    %% Split up conditions
    switch contrast
      case 'congruent'
        idx(1,:) = 11;
        idx(2,:) = 14;
      case 'attended'
        tmpidx{1} = [11 12; 11 13]; % CW:  first row corresponds to left, second to right hemifield
        tmpidx{2} = [13 14; 12 14];
        idx(1,:) = tmpidx{1}(hemi,:);
        idx(2,:) = tmpidx{2}(hemi,:);
    end
    
    for k=1:size(idx,1)
      cfg=[];
      cfg.trials = ismember(data.trialinfo(:,2), idx(k,:));
      dat{k} = ft_selectdata(cfg, data);
    end
    
    % make timelock structure
    cfg=[];
    cfg.keeptrials = 'yes';
    for k=1:numel(dat)
      dat{k} = ft_timelockanalysis(cfg, dat{k});
      ntrl(k) = size(dat{k}.trial,1);
    end
    ntrials = min(ntrl);
    
    %% Noise reduction
    % increase SNR by temporal smoothing
    if do_smooth
      for k=1:numel(dat)
        dat{k}.trial = smoothdata(dat{k}.trial, 3, 'movmean', 4);
      end
    end
      
    % increase SNR by averaging trials randomly
    if do_avgtrials
      groupsize = 5;
      ngroups = floor(ntrials/groupsize);
      for k=1:numel(dat)
        randnum_trlavg = randperm(ntrl(k));
        dat{k} = randavg_trials(dat{k}, ngroups, groupsize, randnum_trlavg);
      end
    else
      ngroups = ntrials;
    end
    
    % prewhiten with covariance matrix
    if do_prewhiten
      dat = prewhiten_data(dat);
    end
    
    
    %% decode
    % loop over time started higher up in case of do_binpertimepoint. if enet
    % is used, don't loop over time (this will be done in parallel jobs).
    if do_enet && ~isempty(tbin)
      allt=tbin;
    elseif do_binpertimepoint
      allt=1;
    else
      allt=1:tlength;
    end
    
    for t=allt
      if do_binpertimepoint
        cnt = itrltime;
      else
        cnt = t;
      end
      
      if do_enet
        % prepare model
        model = dml.enet('family', 'binomial', 'df', 0, 'alpha', 0.1 );
        
        % randomize trials
        for k=1:numel(dat)
          randnum_trlorder(k,:) = randperm(ngroups);
          dat{k}.trial = dat{k}.trial(randnum_trlorder(k,:),:,:);
        end
        
        % loop over trials: leave one trial out decoding
        for itrl=1:ngroups
          [traindata, testdata, traindesign, testdesign] = enet_preparedata(dat, itrl, cnt);
          model = model.train(traindata,traindesign);
          tmpacc = model.test(testdata);
          for k=1:numel(testdesign)
            accuracy(bin, itrl+(k-1)*ngroups,cnt) = tmpacc(k,testdesign(k));
          end
        end
        
      else
        cfgs=[];
        cfgs.mva=  {dml.svm};
        cfgs.method = 'crossvalidate';
        cfgs.statistic = {'accuracy'};
        cfgs.type= 'nfold';
        cfgs.nfolds = ngroups;
        cfgs.design = [ones(ngroups,1); 2*ones(ngroups,1)];
        cfgs.latency = [toi(cnt) toi(cnt+nsample-1)];
        stat(bin,cnt) = ft_timelockstatistics(cfgs,dat{1}, dat{2});
        accuracy(bin,cnt) = stat(bin,cnt).statistic.accuracy;
      end
    end
  end
end

% make filename
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

% save variables
settings = struct('allses', do_allses, 'baselinecorrect', do_baselinecorrect, 'pca', do_pca, 'correcttrials', do_correcttrials, 'phasebin', do_phasebin, 'contrast', contrast, 'avgtrials', do_avgtrials, 'prewhiten', do_prewhiten);
save(filename, 'accuracy','settings')
if do_phasebin
  save(filename, 'phase', 'trl', '-append');
end