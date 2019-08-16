function analysis_decoding(subj, contrast, varargin)
if ~exist('contrast', 'var') || nargin<2,       contrast = 'congruent'; end
if ~exist('subj', 'var')     || nargin<1,       subj=4; end

do_allses          = ft_getopt(varargin, 'do_allses',          true);
do_baselinecorrect = ft_getopt(varargin, 'do_baselinecorrect', false);
do_pca             = ft_getopt(varargin, 'do_pca',             false);
do_correcttrials   = ft_getopt(varargin, 'do_correcttrials',   false);
do_timeresolved    = ft_getopt(varargin, 'do_timeresolved',    false);
do_avgtrials       = ft_getopt(varargin, 'do_avgtrials',       false);
do_prewhiten       = ft_getopt(varargin, 'do_prewhiten',       false);
do_smooth          = ft_getopt(varargin, 'do_smooth',          false);
do_enet            = ft_getopt(varargin, 'do_enet',            false);
do_nestedcv        = ft_getopt(varargin, 'do_nestedcv',        false);
do_phasealign      = ft_getopt(varargin, 'do_phasealign',      false);
bpfreq             = ft_getopt(varargin, 'bpfreq',             [8 12]);
do_phasebin        = ft_getopt(varargin, 'do_phasebin',        false);
do_binpertimepoint = ft_getopt(varargin, 'do_binpertimepoint', false);
randnr             = ft_getopt(varargin, 'randnr',             []);
tbin               = ft_getopt(varargin, 'tbin',               []);


ft_info off

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
  data.grad = tmp{1}.grad;
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
[trl, phase, distance, time] = analysis_alphaphase(data, bpfreq, centerphase);


if do_timeresolved
  toi = -0.1:1/data.fsample:1.2;
  tlength = numel(toi);
  twindow = 1/fs;
  nsample = twindow*fs;
else
  toi = 0.4:1/data.fsample:1.2;
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

if do_phasealign
  alignfreq = mean(bpfreq);
  examplephase = 0:2*pi/(fs/alignfreq):2*pi;
  shorttime = nearest(time, 0.4);
  shortphase = phase(:,shorttime:end);
  shorttime = time(shorttime:end);
  examplephase = repmat(examplephase, [1 ceil(numel(shorttime)/numel(examplephase))+1]);
  shiftframes = floor(fs/alignfreq);
  for k=1:size(phase,1)
    for l=1:shiftframes-1
      tmpphase = examplephase(l+1:l+size(phase,2));
      R(k,l+1) = corr(tmpphase', shortphase(k,:)');
    end
  end
  [~, shiftidx] = max(R');
  shiftidx = shiftidx - 1; % miminal shift is no shift at all.
  for k=1:numel(data.trial)
    data.time{k} = data.time{k} + shiftidx(k)/data.fsample;
  end
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
  trltime = 1:length(time);
else
  trltime=1;
end

%% loop over bins (only >1 when phasebinning)
for bin = 1:size(trl,2)
  % select trial subset (when phasebinning)
  for itrltime = trltime
    cfg=[];
    cfg.trials = trl{itrltime,bin}; 
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
    for k=1:numel(dat)
      trlidx = randperm(ntrl(k));
      dat{k}.trial = dat{k}.trial(trlidx(1:ntrials),:,:);
    end
    
    
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
        dat{k} = randavg_trials(dat{k}, ngroups, groupsize);
      end
    else
      ngroups = ntrials;
    end
    
    % prewhiten with covariance matrix
    if do_prewhiten
      % compute covariance on all data if necessary
      % get covariance matrix based on all trials
      [~, cov] = prewhiten_data(dat);
    else
      cov=[];
    end
    
    %% decode
    % loop over time started higher up in case of do_binpertimepoint. if enet
    % is used, don't loop over time (this will be done in parallel jobs).
    if ~isempty(tbin)
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
        if numel(allt)==1
          cnt=1;
        else
          cnt = t;
        end
      end
      
      % randomize trials
      randnum_trlorder = [];
      for k=1:numel(dat)
        randnum_trlorder(k,:) = randperm(ngroups);
        dat{k}.trial = dat{k}.trial(randnum_trlorder(k,:),:,:);
      end
      
      % prepare model
      if do_enet && do_nestedcv
        accuracy = [];
        nfolds = 5;
        l1_ratio_range = [.1, .3, .5, .7, .9, .95, 1];
        
        groupsize_out = floor(size(dat{1}.trial,1)/nfolds);
        groupsize_out = repmat(groupsize_out, 1, nfolds);
        rem = size(dat{1}.trial,1)-sum(groupsize_out);
        groupsize_out = groupsize_out + [ones(1,rem), zeros(1,nfolds-rem)]; % some folds will have one more trial
        
        for fold_out = 1:nfolds
          for fold_in = 1:nfolds-1
            exclude_trials = sum(groupsize_out(1,1:fold_out))-groupsize_out(fold_out)+1:sum(groupsize_out(1,1:fold_out));
            if do_prewhiten
              tmpcov=cov;
              tmpcov(exclude_trials,:,:)=[];
            end
            tmpdat = dat;
            for ii=1:numel(tmpdat)
              tmpdat{ii}.trial(exclude_trials,:,:) = [];
            end
            groupsize_in = groupsize_out;
            groupsize_in(fold_out) = [];
            [train_in, test_in, traindes_in, testdes_in] = dml_preparedata(tmpdat, sum(groupsize_in(1,1:fold_in))-groupsize_in(fold_in)+1:sum(groupsize_in(1,1:fold_in)), cnt, do_prewhiten, tmpcov);
            for nst = 1:numel(l1_ratio_range)
              model_nstd = dml.enet('family', 'binomial', 'df', 0, 'alpha', l1_ratio_range(nst));
              model_nstd = model_nstd.train(train_in,traindes_in);
              tmp = model_nstd.test(test_in);
              acc_nstd(fold_in, nst) = mean([tmp(1:size(tmp,1)/2,1); tmp(size(tmp,1)/2+1:end,2)]);
            end
          end
          [~, maxidx] = max(mean(acc_nstd,1));
          l1_ratio = l1_ratio_range(maxidx);
          model = dml.enet('family', 'binomial', 'df', 0, 'alpha', l1_ratio);
          
          [train_outer, test_outer, traindes_outer, testdes_outer] = dml_preparedata(dat, sum(groupsize_out(1,1:fold_out))-groupsize_out(fold_out)+1:sum(groupsize_out(1,1:fold_out)), cnt, do_prewhiten, cov);
          model = model.train(train_outer,traindes_outer);
          tmpacc = model.test(test_outer);
          for k=1:size(testdes_outer,1)
            accuracy{bin,cnt}(k+size(testdes_outer,1)*(fold_out-1)) = tmpacc(k,testdes_outer(k));
          end
        end
      else
        if do_enet
          l1_ratio = 0.1;
          model = dml.enet('family', 'binomial', 'df', 0, 'alpha', l1_ratio);
        else
          model = dml.svm;
        end
        
        % loop over trials: leave one trial out decoding
        for itrl=1:ngroups
          itrl
          [traindata, testdata, traindesign, testdesign] = dml_preparedata(dat, itrl, t, do_prewhiten, cov);
          model = model.train(traindata,traindesign);
          tmpacc = model.test(testdata);
          for k=1:numel(testdesign)
            accuracy{bin, cnt}(itrl+(k-1)*ngroups) = tmpacc(k,testdesign(k));
          end
        end
        
        %       else
        %
        %         cfgs=[];
        %         cfgs.mva=  {dml.svm};
        %         cfgs.method = 'crossvalidate';
        %         cfgs.statistic = {'accuracy'};
        %         cfgs.type= 'nfold';
        %         cfgs.nfolds = ngroups;
        %         cfgs.design = [ones(ngroups,1); 2*ones(ngroups,1)];
        %         cfgs.latency = [toi(cnt) toi(cnt+nsample-1)];
        %         stat(bin,cnt) = ft_timelockstatistics(cfgs,dat{1}, dat{2});
        %         accuracy(bin,cnt) = stat(bin,cnt).statistic.accuracy;
        %       end
      end
    end
  end
end
vararg = [];
% make filename
filename = '/project/3011085.02/phasecode/results/';
if do_phasebin
  filename = [filename, 'phasebin_svm/'];
elseif do_enet
  filename = [filename, 'enet/'];
  vararg.l1_ratio = l1_ratio;
else
  filename = [filename, 'svm/'];
end
filename = [filename, sprintf('sub%02d_decoding', subj)];
if do_phasealign
  filename = [filename, '_phasealign'];
end

if exist('tbin', 'var') && ~isempty(tbin)
  filename = [filename, sprintf('_%d', tbin)];
end
if exist('randnr', 'var')
  filename = [filename, sprintf('_%d', randnr)];
end

% save variables
settings = struct('allses', do_allses, 'avgtrials', do_avgtrials, 'baselinecorrect', do_baselinecorrect,'binpertimepoint', do_binpertimepoint,'correcttrials', do_correcttrials, 'enet', do_enet, 'nestedcv', do_nestedcv, 'pca', do_pca, 'phasealign', do_phasealign, 'phasebin', do_phasebin, 'contrast', contrast,'prewhiten', do_prewhiten, 'smooth', do_smooth,'timeresolved', do_timeresolved, 'var', vararg);

save(filename, 'accuracy','settings')
if do_phasebin
  save(filename, 'phase', 'trl', '-append');
end