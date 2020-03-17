datainfo;
addpath([projectdir, 'scripts/gcmi/matlab'])
hemi = 1;

% get classification accuracy and testdata for MEG and eyetracker decoding.
rngseed = rand(1)*10^6;

[accuracy{1},~,~,dat(1)] = analysis_decoding(subj,'attended','hemi', hemi,...
'do_correcttrials', 1, 'do_avgtrials', 1, 'rngseed', rngseed,...
'do_prewhiten', 2, 'f', 4,'nbins', 110, 'groupsize', 5,'nperm', 1,...
'do_randphasebin', 0, 'dosave',0, 'keepfolds', true)

[accuracy{2},~,~,dat(2)] = analysis_decoding(subj,'attended','hemi', hemi,...
'do_correcttrials', 1, 'do_avgtrials', 1, 'rngseed', rngseed,...
'do_prewhiten', 2, 'f', 4,'nbins', 110, 'groupsize', 5,'nperm', 1,...
'do_randphasebin', 0, 'dosave',0, 'keepfolds', true, 'chansel', 'eye')

  for k=1:110
    tmp1 = [];
    tmp2 = [];
    for l=1:5
      tmp1 = [tmp1; dat(1).testdata_fold{1,l,k}];
      x{k} = tmp1;
      
      tmp2 = [tmp2; dat(2).testdata_fold{1,l,k}];
      y{k} = tmp2;
    end
  end
  
  for k=1:110
    I(k,1) = gcmi_cc(x{k}, y{k});
  end
  
  r = corr(I, mean(accuracy{1},2))
  
%%
[accuracy{1},~,~,dat(1)] = analysis_decoding(subj,'attended','hemi', hemi,...
'do_correcttrials', 1, 'do_avgtrials', 1, 'rngseed', rngseed,...
'do_prewhiten', 2, 'f', 4,'nbins', 1, 'groupsize', 5,'nperm', 1,...
'do_randphasebin', 0, 'dosave',0, 'keepfolds', true, 'do_randphasebin', 0);

[~,~,~,dat(2)] = analysis_decoding(subj,'attended','hemi', hemi,...
'do_correcttrials', 1, 'do_avgtrials', 1, 'rngseed', rngseed,...
'do_prewhiten', 2, 'f', 4,'nbins', 1, 'groupsize', 5,'nperm', 1,...
'do_randphasebin', 0, 'dosave',0, 'keepfolds', true, 'chansel', 'eye', 'do_randphasebin', 1);

% concatenate folds and select only trials with correct classiciation in MEG.
correct = cat(1,dat(1).accuracy_fold{:});
incorrect = ~correct;
for k=1:2
  testdata_fold = (cat(1,dat(k).testdata_fold{:}));
  
  data{k} = testdata_fold(find(correct),:);
  data2{k} = testdata_fold(find(incorrect),:);
end

% compute MI
I = gcmi_cc(data{1}, data{2});

%% random permutations
% permute trials in eye data
N = size(data{1},1);
nperm=100;

for k=1:nperm
  idx = randperm(N);
  Iperm(k) = gcmi_cc(data{1}, data{2}(idx,:));
end








