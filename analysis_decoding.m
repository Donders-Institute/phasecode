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
nfolds              = ft_getopt(varargin, 'nfolds',              5);
nperm               = ft_getopt(varargin, 'nperm',               10);
nrandperm           = ft_getopt(varargin, 'nrandperm',           100);

if ~do_randphasebin, nrandperm = 1; end

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
clear tmp

fs = data.fsample;


% Find trialindices per condition
cfg=[];
switch contrast
  case 'congruent'
    idx_ori(1,:) = 11;
    idx_ori(2,:) = 14;
    idx = ismember(data.trialinfo(:,2), idx_ori);
  case 'attended'
    tmpidx{1} = [11 12; 11 13]; % CW:  first row corresponds to left, second to right hemifield
    tmpidx{2} = [13 14; 12 14];
    idx_ori(1,:) = tmpidx{1}(hemi,:);
    idx_ori(2,:) = tmpidx{2}(hemi,:);
    idx = ismember(data.trialinfo(:,2), idx_ori) & data.trialinfo(:,1)==hemi;
  case 'unattended'
    tmpidx{1} = [11 12; 11 13]; % CW:  first row corresponds to left, second to right hemifield
    tmpidx{2} = [13 14; 12 14];
    idx_ori(1,:) = tmpidx{1}(hemi,:);
    idx_ori(2,:) = tmpidx{2}(hemi,:);
    idx = ismember(data.trialinfo(:,2), idx_ori) & data.trialinfo(:,1)~=hemi;
end
cfg.trials = idx;
data = ft_selectdata(cfg, data);

% divide trials into phase bins
centerphase = [0 1/3 2/3 1 4/3 5/3]*pi;%[acos(1), acos(0), acos(-1)];
filename = [projectdir, 'results/phase/', sprintf('sub%02d_phase_%d', subj, f(1))];
load(filename)
phasebin(~cfg.trials,:)=[];
phase(~cfg.trials,:)=[];

nchan = numel(data.label);
ncond = 2;
nbins = numel(unique(phasebin(:)));
accuracy = zeros(nrandperm, nperm, nbins, nfolds);

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


%% Split up conditions
data = data_orig;

for k=1:size(idx_ori,1)
  cfg=[];
  cfg.trials = ismember(data.trialinfo(:,2), idx_ori(k,:));
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
  phasebin_orig{k} = phasebin(cfg.trials,:);
end
clear phasebin

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
end

for irandperm = 1:nrandperm
  for k=1:numel(phasebin_orig)
    if do_randphasebin
      [s1, s2] = size(phasebin_orig{k});
      phasebin{k} = phasebin_orig{k}(randperm(s1),:);
    else
      phasebin{k} = phasebin_orig{k};
    end
    phasebin{k} = reshape(phasebin{k}, [], 1);
  end
  
  %% Select data per phase bin
  for bin=1:nbins
    for k=1:ncond
      bindat{bin,k} = dat{k};
      bindat{bin,k}.time = 1;
      bindat{bin,k}.trial = reshape(permute(bindat{bin,k}.trial,[1 3 2]), [], nchan);
      
      % select samples with particular phase
      sampinbin = phasebin{k}==bin;
      nsamp_cond(bin,k) = sum(sampinbin);
      bindat{bin,k}.trial = bindat{bin,k}.trial(sampinbin,:);
    end
  end
  % NOTE: no sub selection of trials until this point. Data can be used
  % multiple times.
  nsmp = min(nsamp_cond(:));
  bindat_orig = bindat;
  
  % Repeat decoding multiple times with different random averages of trials.
  for iperm=1:nperm
    bindat = bindat_orig;
    % select same amount of data for each condition, and each bin.
    for bin = 1:nbins
      for k=1:2
        tmpidx = randperm(nsamp_cond(bin,k));
        bindat{bin,k}.trial = bindat{bin,k}.trial(tmpidx(1:nsmp),:);
        
        % increase SNR by averaging trials randomly
        groupsize = 10;
        ngroups = floor(nsmp/groupsize);
        bindat{bin,k} = randavg_trials(bindat{bin,k}, ngroups, groupsize);
      end
    end
    
    %% decoding
    % loop over time started higher up in case of do_binpertimepoint. if enet
    % is used, don't loop over time (this will be done in parallel jobs).
    
    % initialize folding parameters
    groupsize_fold = floor(ngroups/nfolds);
    groupsize_fold = repmat(groupsize_fold, 1, nfolds);
    rem = ngroups-sum(groupsize_fold);
    groupsize_fold = groupsize_fold + [ones(1,rem), zeros(1,nfolds-rem)];
    
    for bin = 1:nbins
      % initialize model
      model = dml.svm;
      
      % n-fold cross validation
      for ifold=1:nfolds
        % select data per fold and pre-whiten
        [traindata, testdata, traindesign, testdesign] = dml_preparedata(bindat(bin,:), sum(groupsize_fold(1,1:ifold))-groupsize_fold(ifold)+1:sum(groupsize_fold(1,1:ifold)), 1, do_prewhiten);
        
        model = model.train(traindata,traindesign);
        primal{iperm, bin,ifold} = model.primal;
        tmpacc = model.test(testdata);
        for k=1:numel(testdesign)
          tmpacc2(k) = tmpacc(k,testdesign(k));
        end
        accuracy(irandperm, iperm, bin, ifold) = mean(tmpacc2);
        clear tmpacc2
      end
    end
    clear bindat
  end
end
if size(accuracy,1)==1, accuracy = squeeze(accuracy); end

%% save
vararg = [];
% make filename
filename = [projectdir 'results/collapsephasebin/'];

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

if exist('randnr', 'var') && ~isempty(randnr)
  filename = [filename, sprintf('_%d', randnr)];
end

% save variables
settings = struct('avgtrials', do_avgtrials, 'correcttrials', do_correcttrials, 'contrast', contrast,'prewhiten', do_prewhiten, 'var', vararg, 'nperm',nperm, 'nrandperm', nrandperm);
if ~exist('primal', 'var'), primal=[]; end

save(filename, 'accuracy','settings', 'primal')
