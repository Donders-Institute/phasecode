function [accuracy,centerphase,primal] = analysis_decoding(subj, contrast, varargin)
if ~exist('contrast', 'var') || nargin<2,       contrast = 'congruent'; end
if ~exist('subj', 'var')     || nargin<1,       subj=4; end

do_correcttrials    = ft_getopt(varargin, 'do_correcttrials',    false);
do_avgtrials        = ft_getopt(varargin, 'do_avgtrials',        false);
do_prewhiten        = ft_getopt(varargin, 'do_prewhiten',        0);
randnr              = ft_getopt(varargin, 'randnr',              []);
rngseed             = ft_getopt(varargin, 'rngseed',             []);
hemi                = ft_getopt(varargin, 'hemi',                1);
f                   = ft_getopt(varargin, 'f',                   10);
do_randphasebin     = ft_getopt(varargin, 'do_randphasebin',     false);
nfolds              = ft_getopt(varargin, 'nfolds',              5);
nperm               = ft_getopt(varargin, 'nperm',               10);
nrandperm           = ft_getopt(varargin, 'nrandperm',           100);
nbins               = ft_getopt(varargin, 'nbins',               12);
groupsize           = ft_getopt(varargin, 'groupsize',           10);
binoverlap          = ft_getopt(varargin, 'binoverlap',          false);
tmpfilename         = ft_getopt(varargin, 'tmpfilename',         []);
dosave              = ft_getopt(varargin, 'dosave',              true);

if ~do_randphasebin, nrandperm = 1; end
if isempty(rngseed)
  rng('shuffle'); rngseed = rand(1)*10^6;
end
rng(rngseed, 'twister');

ft_info off

% load data
datainfo;
cfg=[];
cnt=1;
for ses=subjects(subj).validsessions
  filename = [datadir, sprintf('sub%02d/meg%02d/sub%02d-meg%02d_cleandata.mat', subj, ses, subj, ses)];
  tmp{cnt} = load(filename, 'data');
  tmp{cnt} = removefields(tmp{cnt}.data, 'elec');
  cnt=cnt+1;
end
data = ft_appenddata([], tmp{:});
clear tmp

fs = data.fsample;


% Find trialindices per condition
[data, idx, idx_ori] = phasecode_select_condition(data, contrast, hemi);

% divide trials into phase bins
filename = [projectdir, 'results/phase/', sprintf('sub%02d_phase_%d', subj, f(1))];
load(filename)
[phasebin, phase, ~, time]  = analysis_phase(subj, f(1), [], 3, false);
% if strcmp(contrast, 'unattended')
%   phasebin = phasebin{hemi};
%   phase = phase{hemi};
% else
  phasebin = phasebin{mod(hemi,2)+1}; % for attention left, take right hemisphere phase.
  phase = phase{mod(hemi,2)+1};


% end
phasebin(~idx,:)=[];
phase(~idx,:)=[];

if 1 % only use data after the main ERF (post 400 ms)
  time(1:80)=[];
  phasebin(:,1:80)=[];
  phase(:,1:80) = [];
end 

nchan = numel(data.label);
ncond = 2;
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
  phase_orig{k}    = phase(cfg.trials,:);
end
clear phasebin phase

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
      % circularly shift each row by a random amount, maximal 1 period
      % forward or backward in time.
      [s1, s2] = size(phasebin_orig{k});
      shufvec     = randi([-floor(fs/f) floor(fs/f)],s1,1);
      phasebin{k} = circshift_columns(transpose(phasebin_orig{k}), shufvec)';
      phase{k}    = circshift_columns(transpose(phase_orig{k}), shufvec)';
    else
      phasebin{k} = phasebin_orig{k};
      phase{k}    = phase_orig{k};
    end
    phasebin{k} = reshape(phasebin{k}, [], 1);
    phase{k}    = reshape(phase{k}, [] ,1);
  end
  
  %% Select data per phase bin
  
  if binoverlap
    centerphase = (0:1/nbins:(nbins-1)/nbins)*pi;
    for k=1:numel(phase)
      for l=1:numel(phase{k})
        d{k}(l,:) = abs(circ_dist(phase{k}(l), centerphase));
      end
      tmpn = size(d{k},1);
      tmpn_perbin = tmpn/5;
      usepercentage = 0.6;
      usensmp = floor(tmpn_perbin*usepercentage);
      [~, lowestd{k}] = sort(d{k}, 1, 'ascend');
      usesmp{k} = lowestd{k}(1:usensmp,:);
    end
  else
    centerphase = linspace(0,2*pi,nbins+1);
    centerphase = [centerphase(end-1)-2*pi centerphase];
  end
  for bin=1:nbins
    for k=1:ncond
      bindat{bin,k} = dat{k};
      bindat{bin,k}.time = 1;
      bindat{bin,k}.trial = reshape(permute(bindat{bin,k}.trial,[1 3 2]), [], nchan);
      if binoverlap
        bindat{bin,k}.trial = bindat{bin,k}.trial(usesmp{k}(:,bin),:);
        nsamp_cond(bin,k) = numel(usesmp{k}(:, bin));
      else
      % select samples with particular phase
      if bin==1
        % deal with wrap around 0
        sampinbin = phase{k}-2*pi>=centerphase(bin)|phase{k}<centerphase(bin+2);
      else
        sampinbin = phase{k}>=centerphase(bin)&phase{k}<centerphase(bin+2);
      end
      
      nsamp_cond(bin,k) = sum(sampinbin);
      bindat{bin,k}.trial = bindat{bin,k}.trial(sampinbin,:);
      end
    end
  end
  % NOTE: no sub selection of trials until this point. Data can be used
  % multiple times.
  nsmp = min(nsamp_cond(:));
  bindat_orig = bindat;
  
  % Repeat decoding multiple times with different random averages of trials.
  ngroups   = floor(nsmp/groupsize);
  x = repmat(1:ngroups,[groupsize,1]); x = x(:);
  y = 1:(groupsize*ngroups); y = y(:);
  z = ones(size(y))./groupsize;
  for iperm=1:nperm
    fprintf('computing permutation %d/%d\n',iperm,nperm);
    bindat = bindat_orig;
    % select same amount of data for each condition, and each bin.

%{
%% the commented out section is on average 5 times as slow as the sparse multiplication
%
%     % select same amount of data for each condition, and each bin.
%     for bin = 1:nbins
%       for k=1:2
%         tmpidx = randperm(nsamp_cond(bin,k));
%         bindat{bin,k}.trial = bindat{bin,k}.trial(tmpidx(1:nsmp),:);
%         
%         % increase SNR by averaging trials randomly
%         groupsize = 10;
%         ngroups = floor(nsmp/groupsize);
%         bindat{bin,k} = randavg_trials(bindat{bin,k}, ngroups, groupsize);
%       end
%     end
%}
    for bin = 1:nbins
      P    = sparse(x,y,z,ngroups,size(bindat{bin,1}.trial,1));
      tmp1 = P(:,randperm(size(P,2)))*bindat{bin,1}.trial;
      P    = sparse(x,y,z,ngroups,size(bindat{bin,2}.trial,1));
      tmp2 = P(:,randperm(size(P,2)))*bindat{bin,2}.trial;
       
      if do_prewhiten==2
        tmp     = [tmp1;tmp2];
        tmp     = tmp-mean(tmp);
        [u,s,v] = svd(tmp'*tmp);
        diagS   = diag(s);
        thr     = 0.99;
        
        sel     = find(cumsum(diagS)./sum(diagS)<=thr);
        P       = diag(1./sqrt(diagS(sel)))*u(:,sel)';
        tmp1    = tmp1*P';
        tmp2    = tmp2*P';
      end
      bindat{bin,1}.trial = tmp1;
      bindat{bin,2}.trial = tmp2;  
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
      % n-fold cross validation
      for ifold=1:nfolds
        % select data per fold and pre-whiten
        indx_testtrials = sum(groupsize_fold(1,1:ifold))-groupsize_fold(ifold)+1:sum(groupsize_fold(1,1:ifold));
        [traindata, testdata, traindesign, testdesign, tmpP] = dml_preparedata(bindat(bin,:), indx_testtrials, 1, do_prewhiten);
        if ~isempty(tmpP), P=tmpP; end
        
        % initialize model
        model          = dml.svm;
        model.kernel   = 'precomp';
        model.issqrtmK = true;
        model.distance = true;
        if ifold == 1
          Cparam = 0.1.*mean(sum([traindata;testdata].^2,2));
        end
        model.C      = Cparam;
        
        [u,s,v]      = svd(traindata,'econ');
        model.Ktrain = u*s;
        
        model = model.train(traindata,traindesign);
        primal{iperm, bin,ifold} = model.primal;
        tmpacc = model.test(testdata);
        accuracy(irandperm, iperm, bin, ifold) = mean(double(double(tmpacc(:,1)>0)==mod(testdesign,2))); % this line works in a 2 class classification problem
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
filename = [filename, sprintf('sub%02d/phase3/sub%02d_decoding', subj, subj)];
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
settings = struct('avgtrials', do_avgtrials, 'correcttrials', do_correcttrials, 'contrast', contrast,'prewhiten', do_prewhiten>0, 'var', vararg, 'nperm',nperm, 'nrandperm', nrandperm);
if ~exist('primal', 'var') || do_randphasebin, primal=[]; end
rnumb=rng;
if dosave
    if isempty(tmpfilename)
      save(filename, 'accuracy','settings', 'primal', 'P', 'rnumb')
    else
      save(tmpfilename, 'accuracy', 'rnumb')
    end
end