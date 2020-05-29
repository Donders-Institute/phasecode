function control_analysis_eye(subj, hemi)
datainfo;
rng('default')
rngseed = 310520;

%% eye
% get data
[~,~,~,~,returndata_eye] = analysis_decoding(subj,'attended','hemi', hemi,...
  'do_correcttrials', 1, 'do_avgtrials', 1, 'rngseed', rngseed,...
  'do_prewhiten', 2, 'f', 4,'nbins', 1, 'groupsize', 5, 'chansel', 'eye',...
  'do_randphasebin',0, 'randnr', [], 'nrandperm',1, 'nperm', 100, 'dosave', 0, 'do_controleye',1);
nfolds=10;
ngroups = returndata_eye.ngroups;
groupsize_fold = floor(ngroups/nfolds);
groupsize_fold = repmat(groupsize_fold, 1, nfolds);
rem = ngroups-sum(groupsize_fold);
groupsize_fold = groupsize_fold + [ones(1,rem), zeros(1,nfolds-rem)];

for ifold=1:nfolds
  % select data per fold and pre-whiten
  indx_testtrials = sum(groupsize_fold(1,1:ifold))-groupsize_fold(ifold)+1:sum(groupsize_fold(1,1:ifold));
  [traindata{ifold}, testdata{ifold}, traindesign{ifold}, testdesign{ifold}] = dml_preparedata(returndata_eye.bindat, indx_testtrials, 1, 2);
  trial_index_orig{ifold} = 1:size(traindata{ifold},1);
end
% trial_index = trial_index_orig;

idx_all_trials = 1:ngroups;
idx_traindata_trials = [];

data_eye = returndata_eye.bindat;


T = 1;
N=0;
while T==1
  acc=[];
  acc0=[];
  clear dist which_correct
  if N==0
    ngroups = returndata_eye.ngroups;
    for k=1:2
      idx_leftover_trials{k} = idx_all_trials;
      idx_rem_trials{k} = [];
    end
  else
    ngroups = min(size(data_eye{1}.trial,1),size(data_eye{2}.trial,1));
    for k=1:2
      if size(data_eye{k}.trial,1)>ngroups
        data_eye{k}.trial(ngroups+1:end,:) = [];
        idx_rem_trials{k} = [idx_rem_trials{k}, idx_leftover_trials{k}(ngroups+1:end)];
        idx_leftover_trials{k}(ngroups+1:end) = [];
      end
    end
  end

groupsize_fold = floor(ngroups/nfolds);
groupsize_fold = repmat(groupsize_fold, 1, nfolds);
rem = ngroups-sum(groupsize_fold);
groupsize_fold = groupsize_fold + [ones(1,rem), zeros(1,nfolds-rem)];
randomizer = randperm(ngroups);
for ifold=1:nfolds
  % select data per fold and pre-whiten
  indx_testtrials = randomizer(sum(groupsize_fold(1:ifold-1))+1:sum(groupsize_fold(1:ifold)));
  idx_traindata_trials{ifold} = setdiff(1:ngroups, indx_testtrials);
  [traindata{ifold}, testdata{ifold}, traindesign{ifold}, testdesign{ifold}] = dml_preparedata(data_eye, indx_testtrials, 1, 2);
  trial_index_orig{ifold} = 1:size(traindata{ifold},1);
end
  
  for ifold=1:nfolds
    % initialize model
    model          = dml.svm;
    model.kernel   = 'precomp';
    model.issqrtmK = true;
    model.distance = true;
    if ifold == 1
      Cparam = 0.1.*mean(sum([traindata{ifold};testdata{ifold}].^2,2));
    end
    model.C      = Cparam;
    [u,s,v]      = svd(traindata{ifold},'econ');
    model.Ktrain = u*s;
    
    
    model = model.train(traindata{ifold},traindesign{ifold});
    model0 = model;
    model0 = model0.train(traindata{ifold},traindesign{ifold}(randperm(numel(traindesign{ifold}))));
    dual{ifold} = model.dual(1:end-1);
    dist{ifold} = model.test(testdata{ifold});
    which_correct{ifold} = double(double(dist{ifold}(:,1)>0)==mod(testdesign{ifold},2));
    acc(ifold) = mean(which_correct{ifold}); %
    tmpacc0{ifold} = model0.test(testdata{ifold});
    acc0(ifold) = mean(double(double(tmpacc0{ifold}(:,1)>0)==mod(testdesign{ifold},2))); %
  end
  T=ttest(acc,acc0);
  if T==1
    if N==0
      acc_eye_pre = acc;
      acc0_eye_pre = acc0;
      [T_eye_pre.H, T_eye_pre.P, T_eye_pre.CI, T_eye_pre.STATS] = ttest(acc_eye_pre,acc0_eye_pre);
    end
    N=N+1;
    for ifold=1:nfolds
      n=size(trial_index{ifold},2);
      dual1 = dual{ifold}(1:n/2,1);
      dual2 = dual{ifold}(n/2+1:end,1);
      % find the most discriminating observation for each class
      [~, ix1] = max(abs(dual1));
      [~, ix2] = max(abs(dual2));
      
      tmpidx1(ifold) = idx_traindata_trials{ifold}(ix1); % this is the trial index from eye_data
      tmpidx2(ifold) = idx_traindata_trials{ifold}(ix2);
    end
    % remove the trials with the largest dual weights from all folds
    idx_rem_trials{1} = [idx_rem_trials{1}, idx_leftover_trials{1}(unique(tmpidx1))]; % remember which original pseudo trial we throw away
    data_eye{1}.trial(unique(tmpidx1),:) = []; % remove from data
    idx_leftover_trials{1}(unique(tmpidx1)) = []; % remove from trial index pool
    idx_rem_trials{2} = [idx_rem_trials{2}, idx_leftover_trials{2}(unique(tmpidx2))];
    data_eye{2}.trial(unique(tmpidx2),:) = [];
    idx_leftover_trials{2}(unique(tmpidx2)) = [];
  else
    acc_eye = acc;
    acc0_eye = acc0;
    [T_eye.H, T_eye.P, T_eye.CI, T_eye.STATS] = ttest(acc_eye,acc0_eye);
  end
end

percentage_removed = 100-ngroups./returndata_eye.ngroups*100;
n_trials_removed = returndata_eye.ngroups-ngroups;


%% MEG
[~,~,~,~,returndata_MEG] = analysis_decoding(subj,'attended','hemi', hemi,...
  'do_correcttrials', 1, 'do_avgtrials', 1, 'rngseed', rngseed,...
  'do_prewhiten', 2, 'f', 4,'nbins', 1, 'groupsize', 5, 'chansel', 'MEG',...
  'do_randphasebin',0, 'randnr', [], 'nrandperm',1, 'nperm', 100, 'dosave', 0, 'do_controleye',1);

data_meg = returndata_MEG.bindat;
for ii = 1:2 % do it once for all data, once without the trials discarded from eye decoding
  
  if ii==2
    for k=1:2 % for both classes
      data_meg{k}.trial(idx_rem_trials{k},:)=[];
    end
  end
  
  randomizer = randperm(size(data_meg{1}.trial,1));
  for ifold=1:nfolds
    % select data per fold and pre-whiten
    indx_testtrials = randomizer(sum(groupsize_fold(1:ifold-1))+1:sum(groupsize_fold(1:ifold)));
    [traindata{ifold}, testdata{ifold}, traindesign{ifold}, testdesign{ifold}] = dml_preparedata(data_meg, indx_testtrials, 1, 2);
  end
  
  
  for ifold=1:nfolds
    % initialize model
    model          = dml.svm;
    model.kernel   = 'precomp';
    model.issqrtmK = true;
    model.distance = true;
    if ifold == 1
      Cparam = 0.1.*mean(sum([traindata{ifold};testdata{ifold}].^2,2));
    end
    model.C      = Cparam;
    [u,s,v]      = svd(traindata{ifold},'econ');
    model.Ktrain = u*s;
    
    
    model = model.train(traindata{ifold},traindesign{ifold});
    model0 = model;
    model0 = model0.train(traindata{ifold},traindesign{ifold}(randperm(numel(traindesign{ifold}))));
    dual{ifold} = model.dual(1:end-1);
    dist{ifold} = model.test(testdata{ifold});
    which_correct{ifold} = double(double(dist{ifold}(:,1)>0)==mod(testdesign{ifold},2));
    acc(ifold) = mean(which_correct{ifold}); %
    tmpacc0{ifold} = model0.test(testdata{ifold});
    acc0(ifold) = mean(double(double(tmpacc0{ifold}(:,1)>0)==mod(testdesign{ifold},2))); %
  end
  clear T
  if ii==1
    acc_MEG_pre = acc;
    acc0_MEG_pre = acc0;
    [T_MEG_pre.H, T_MEG_pre.P, T_MEG_pre.CI, T_MEG_pre.STATS] = ttest(acc_MEG_pre,acc0_MEG_pre);
  else
    acc_MEG = acc;
    acc0_MEG = acc0;
    [T_MEG.H, T_MEG.P, T_MEG.CI, T_MEG.STATS] = ttest(acc_MEG,acc0_MEG);
  end
end

if ~exist('T_eye_pre', 'var'), T_eye_pre = []; end
if ~exist('acc_eye_pre', 'var'), acc_eye_pre = []; end
if ~exist('acc0_eye_pre', 'var'), acc0_eye_pre = []; end

save([results_dir, 'control/', sprintf('sub%02d_analysis_eye_hemi%d', subj, hemi)], 'N','T_eye_pre','T_eye', 'acc0_eye_pre', 'acc_eye_pre', 'acc0_eye', 'acc_eye', 'T_MEG', 'acc_MEG', 'acc0_MEG','T_MEG_pre', 'acc_MEG_pre', 'acc0_MEG_pre','percentage_removed', 'n_trials_removed')
