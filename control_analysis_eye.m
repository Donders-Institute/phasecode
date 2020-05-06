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
trial_index = trial_index_orig;

T = 1;
N=0;
while T==1
  acc=[];
  acc0=[];
  clear dist which_correct
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
    N=N+1;
    for ifold=1:nfolds
      n=size(trial_index{ifold},2);
      dual1 = dual{ifold}(1:n/2,1);
      dual2 = dual{ifold}(n/2+1:end,1);
      % find the most discriminating observation for each class
      [~, ix1] = max(abs(dual1));
      [~, ix2] = max(abs(dual2));
      ix2 = ix2+n/2; % add up the amount of trials from class 1
      % remove these trials from the traindata, traindesign, and trial index
      traindata{ifold}([ix1, ix2],:)=[];
      traindesign{ifold}([ix1 ix2],:)=[];
      trial_index{ifold}(:,[ix1 ix2])=[];
    end
  else
    acc_eye = acc;
    acc0_eye = acc0;
    [T_eye.H, T_eye.P, T_eye.CI, T_eye.STATS] = ttest(acc_eye,acc0_eye);
  end
end

percentage_removed = 100-numel(trial_index{1})./numel(trial_index_orig{1})*100;
for ifold=1:nfolds
  remove_trials{ifold} = setdiff(trial_index_orig{ifold}, trial_index{ifold});
end


%% MEG
[~,~,~,~,returndata_MEG] = analysis_decoding(subj,'attended','hemi', hemi,...
  'do_correcttrials', 1, 'do_avgtrials', 1, 'rngseed', rngseed,...
  'do_prewhiten', 2, 'f', 4,'nbins', 1, 'groupsize', 5, 'chansel', 'MEG',...
  'do_randphasebin',0, 'randnr', [], 'nrandperm',1, 'nperm', 100, 'dosave', 0, 'do_controleye',1);

for ifold=1:nfolds
  % select data per fold and pre-whiten
  indx_testtrials = sum(groupsize_fold(1,1:ifold))-groupsize_fold(ifold)+1:sum(groupsize_fold(1,1:ifold));
  [traindata{ifold}, testdata{ifold}, traindesign{ifold}, testdesign{ifold}] = dml_preparedata(returndata_MEG.bindat, indx_testtrials, 1, 2);
  traindata{ifold}(remove_trials{ifold},:)=[];
  traindesign{ifold}(remove_trials{ifold},:)=[];
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
acc_MEG = acc;
acc0_MEG = acc0;
[T_MEG.H, T_MEG.P, T_MEG.CI, T_MEG.STATS] = ttest(acc_MEG,acc0_MEG);

save([results_dir, 'control/', sprintf('sub%02d_analysis_eye_hemi%d', subj, hemi)], 'T_eye', 'acc0_eye', 'acc_eye', 'T_MEG', 'acc_MEG', 'acc0_MEG','percentage_removed')
