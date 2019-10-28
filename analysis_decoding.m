function analysis_decoding(subj, contrast, varargin)
if ~exist('contrast', 'var') || nargin<2,       contrast = 'congruent'; end
if ~exist('subj', 'var')     || nargin<1,       subj=4; end

do_correcttrials    = ft_getopt(varargin, 'do_correcttrials',    false);
do_avgtrials        = ft_getopt(varargin, 'do_avgtrials',        false);
do_prewhiten        = ft_getopt(varargin, 'do_prewhiten',        false);
randnr              = ft_getopt(varargin, 'randnr',              []);
hemi                = ft_getopt(varargin, 'hemi',                1);
f                   = ft_getopt(varargin, 'f',                   10);
do_randphasebin     = ft_getopt(varargin, 'do_randphasebin',     false);

if do_collapsephasebin, do_timeresolved = false; end

ft_info off

rng('shuffle');

% load data
datainfo;
  cnt=1;
  for ses=subjects(subj).validsessions
    filename = [datadir, sprintf('sub%02d/meg%02d/sub%02d-meg%02d_cleandata.mat', subj, ses, subj, ses)];
    tmp{cnt} = load(filename, 'data');
    tmp{cnt} = tmp{cnt}.data;
    cnt=cnt+1;
  end
  data = ft_appenddata([], tmp{:});
  data.grad = tmp{1}.grad;
  clear tmp

fs = data.fsample;

% Split up conditions
cfg=[];
switch contrast
  case 'congruent'
    idx(1,:) = 11;
    idx(2,:) = 14;
    cfg.trials = ismember(data.trialinfo(:,2), idx);
  case 'attended'
    tmpidx{1} = [11 12; 11 13]; % CW:  first row corresponds to left, second to right hemifield
    tmpidx{2} = [13 14; 12 14];
    idx(1,:) = tmpidx{1}(hemi,:);
    idx(2,:) = tmpidx{2}(hemi,:);
    cfg.trials = ismember(data.trialinfo(:,2), idx) & data.trialinfo(:,1)==hemi;
  case 'unattended'
    tmpidx{1} = [11 12; 11 13]; % CW:  first row corresponds to left, second to right hemifield
    tmpidx{2} = [13 14; 12 14];
    idx(1,:) = tmpidx{1}(hemi,:);
    idx(2,:) = tmpidx{2}(hemi,:);
    cfg.trials = ismember(data.trialinfo(:,2), idx) & data.trialinfo(:,1)~=hemi;
end
data = ft_selectdata(cfg, data);

  % divide trials into phase bins
  centerphase = [0 1/3 2/3 1 4/3 5/3]*pi;%[acos(1), acos(0), acos(-1)];
  filename = [projectdir, 'results/phase/', sprintf('sub%02d_phase_%d', subj, f(1))];
  load(filename)
  phasebin(~cfg.trials,:)=[];
  phase(~cfg.trials,:)=[];
  if do_randphasebin
    phasebin = reshape(phasebin(randperm(size(phasebin,1)),:), [size(phasebin,1) size(phasebin,2)]);
  end

% select time window of interest.
cfg=[];
cfg.toilim = [time(1) time(end)];
data = ft_redefinetrial(cfg, data);

% only select validly cued and correct trials
if do_correcttrials
  cfg=[];
  cfg.trials = (data.trialinfo(:,7)==1) & (data.trialinfo(:,1)==data.trialinfo(:,4));
  data = ft_selectdata(cfg, data);
  rmtrials = find(~cfg.trials);
  phase(rmtrials,:,:)=[];
  phasebin(rmtrials,:,:)=[];
end
data_orig=data;
nbins = numel(unique(phasebin(:)));

cov=[];

%% Prepare all data.
% loop over bins (only >1 when phasebinning)
for bin = 1:nbins
  % select trial subset (when phasebinning)
      data = data_orig;
         
    %% Split up conditions
    for k=1:size(idx,1)
      cfg=[];
      cfg.trials = ismember(data.trialinfo(:,2), idx(k,:));
      if iscell(data.trial)
        % use the following code to speed op computation
        dat{k} = data;
        dat{k}.trial = dat{k}.trial(cfg.trials);
        dat{k}.time = dat{k}.time(cfg.trials);
        dat{k}.trialinfo = dat{k}.trialinfo(cfg.trials,:);
        dat{k}.sampleinfo = dat{k}.sampleinfo(cfg.trials,:);
      else
        dat{k} = ft_selectdata(cfg, data);
      end
      phasebin_cond{k} = phasebin(cfg.trials,:);
    end
    
    % make timelock structure
    cfg=[];
    cfg.keeptrials = 'yes';
    for k=1:numel(dat)
      try % save some time by not using fieldtrip function (only works if all trials have the same time axis)
        dat{k}.trial = permute(cat(3, dat{k}.trial{:}), [3 1 2]);
        dat{k}.time = dat{k}.time{1};
        dat{k}.dimord = 'rpt_chan_time';
      catch
        dat{k} = ft_timelockanalysis(cfg, dat{k});
      end
      ntrl(k) = size(dat{k}.trial,1);
    end
    ntrials = min(ntrl);
    for k=1:numel(dat)
      trlidx = randperm(ntrl(k));
      dat{k}.trial = dat{k}.trial(trlidx(1:ntrials),:,:);
      phasebin_cond{k} = phasebin_cond{k}(trlidx(1:ntrials),:);
    end
    
    
    %% Noise reduction
    
      nchan = size(dat{1}.trial,2);
      for k=1:numel(dat)
        dat{k}.time = 1;
        dat{k}.trial = reshape(permute(dat{k}.trial,[1 3 2]), [], nchan);
        phasebin_cond{k} = reshape(phasebin_cond{k}, [],1);
      end
      for k=1:numel(phasebin_cond)
        nsamp_cond(k) = sum(phasebin_cond{k}==bin);
      end
      nsamp = min(nsamp_cond);
      
      for k=1:numel(phasebin_cond)
        sampinbin = phasebin_cond{k}==bin; % select samples with particular phase
        dat{k}.trial = dat{k}.trial(sampinbin,:); 
        phasebin_cond{k} = phasebin_cond{k}(sampinbin,:); % do the same with sample matrix
        r = randperm(nsamp_cond(k));
        dat{k}.trial = dat{k}.trial(r(1:nsamp),:); % select the same amount of samples for all conditions
        phasebin_cond{k} = phasebin_cond{k}(r(1:nsamp));
      end
      ntrials = numel(phasebin_cond{1});
    
      
      % FIXME no sub selection of trials until this point. Data can be
      % used multiple times.
      
    % increase SNR by averaging trials randomly
    if do_avgtrials
        groupsize = 10;
      ngroups = floor(ntrials/groupsize);
      for k=1:numel(dat)
        dat{k} = randavg_trials(dat{k}, ngroups, groupsize);
      end
    else
      ngroups = ntrials;
    end
    
    loopdata(1, bin).ngroups = ngroups;
    loopdata(1, bin).dat = dat;
    loopdata(1, bin).ntrials = ntrials;
    loopdata(1, bin).phasebin_cond = phasebin_cond;
    loopdata(1, bin).nsamp_cond = nsamp_cond;
    
    % prewhiten with covariance matrix
    % this is moved to dml_preparedata because covariance should be
    % computed over trials, not time. Thus, covariance should be computed
    % for every fold (Guggenmos et al 2018).    
end

%% select same amount of data for each bin
  for bin = 1:nbins
    ngroups_all(bin) = loopdata(1, bin).ngroups;  
  end
  ngroups_all = min(ngroups_all);
  
for bin = 1:nbins
  for k=1:numel(loopdata(1, bin).dat)
    idx = randperm(loopdata(1, bin).ngroups);
    loopdata(1, bin).dat{k}.trial = loopdata(1, bin).dat{k}.trial(idx(1:ngroups_all),:);
  end
end 

%% decode
% loop over time started higher up in case of do_binpertimepoint. if enet
% is used, don't loop over time (this will be done in parallel jobs).
for bin = 1:nbins
  % select trial subset (when phasebinning)
    
    clear dat;
    dat = loopdata(1, bin).dat;
    
    if ~isempty(tbin)
      allt=tbin;
    elseif do_binpertimepoint || do_collapsephasebin
      allt=1;
    else
      allt=1:tlength;
    end
    
     
      % prepare model
          nfolds = 5;
          groupsize = floor(size(dat{1}.trial,1)/nfolds);
          groupsize = repmat(groupsize, 1, nfolds);
          rem = size(dat{1}.trial,1)-sum(groupsize);
          groupsize = groupsize + [ones(1,rem), zeros(1,nfolds-rem)];
          
          model = dml.svm;
        
        % loop over trials: leave one trial out decoding
        for itrl=1:nfolds
            [traindata, testdata, traindesign, testdesign] = dml_preparedata(dat, sum(groupsize(1,1:itrl))-groupsize(itrl)+1:sum(groupsize(1,1:itrl)), cnt, do_prewhiten);
          
            model = model.train(traindata,traindesign);
          primal{bin}(itrl, cnt, :) = model.primal;
          tmpacc = model.test(testdata);
          for k=1:numel(testdesign)
              tmpacc2(k) = tmpacc(k,testdesign(k));
          end
            accuracy(bin, itrl) = mean(tmpacc2);
            clear tmpacc2
        end
end
vararg = [];
% make filename
filename = [project_dir 'results/collapsephasebin/'];

switch contrast
  case 'congruent'
    filename = [filename, 'congruent/'];
  case 'attended'
    filename = [filename, 'attended/'];
  case 'unattended'
    filename = [filename, 'unattended/'];
end
filename = [filename, sprintf('sub%02d_decoding', subj)];
if do_randphasebin
  filename = [filename, '_rand'];
end
if strcmp(contrast, 'attended') || strcmp(contrast, 'unattended')
  filename = [filename, sprintf('_hemi%d', hemi)];
end
  filename = [filename, sprintf('_f%d', f)];

if exist('randnr', 'var')
  filename = [filename, sprintf('_%d', randnr)];
end

% save variables
settings = struct('avgtrials', do_avgtrials, 'correcttrials', do_correcttrials, 'contrast', contrast,'prewhiten', do_prewhiten, 'var', vararg);
if ~exist('primal', 'var'), primal=[]; end

save(filename, 'accuracy','settings', 'primal')
