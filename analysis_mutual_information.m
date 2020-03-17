function [r, r_rand, MI, MI_rand] = analysis_mutual_information(subj)
datainfo;
addpath([projectdir, 'scripts/gcmi/matlab'])
hemis = 1:2;

% get classification accuracy and testdata for MEG and eyetracker decoding.
rngseed = rand(1)*10^6;
for h=hemis
  [accuracy{1,h},~,~,dat(1,h)] = analysis_decoding(subj,'attended','hemi', h,...
    'do_correcttrials', 1, 'do_avgtrials', 1, 'rngseed', rngseed,...
    'do_prewhiten', 2, 'f', 4,'nbins', 110, 'groupsize', 5,'nperm', 1,...
    'do_randphasebin', 0, 'dosave',0, 'keepfolds', true, 'timedecoding', true);
  
  [accuracy{2,h},~,~,dat(2,h)] = analysis_decoding(subj,'attended','hemi', h,...
    'do_correcttrials', 1, 'do_avgtrials', 1, 'rngseed', rngseed,...
    'do_prewhiten', 2, 'f', 4,'nbins', 110, 'groupsize', 5,'nperm', 1,...
    'do_randphasebin', 0, 'dosave',0, 'keepfolds', true, 'chansel', 'eye','timedecoding', true);
  
  for k=1:110
    tmp1 = [];
    tmp2 = [];
    for l=1:5
      tmp1 = [tmp1; dat(1,h).testdata_fold{1,l,k}];
      x{k} = tmp1;
      
      tmp2 = [tmp2; dat(2,h).testdata_fold{1,l,k}];
      y{k} = tmp2;
    end
  end
  
  for k=1:110
    MI(k,h) = gcmi_cc(x{k}, y{k});
  end
  
  r(h) = corr(MI(:,h), mean(accuracy{1,h},2));
  
  %% random permutations
  [accrand{1,h},~,~,datrand{1,h}] = analysis_decoding(subj,'attended','hemi', h,...
    'do_correcttrials', 1, 'do_avgtrials', 1, 'rngseed', rngseed,...
    'do_prewhiten', 2, 'f', 4,'nbins', 110, 'groupsize', 5,'nperm', 1,...
    'do_randphasebin', 1, 'dosave',0, 'keepfolds', true, 'chansel', 'MEG','timedecoding',1, 'nrandperm',100);
  accrand{1,h} = squeeze(mean(accrand{1,h},4));
  
  [accrand{2,h},~,~,datrand{2,h}] = analysis_decoding(subj,'attended','hemi', h,...
    'do_correcttrials', 1, 'do_avgtrials', 1, 'rngseed', rngseed,...
    'do_prewhiten', 2, 'f', 4,'nbins', 110, 'groupsize', 5,'nperm', 1,...
    'do_randphasebin', 1, 'dosave',0, 'keepfolds', true, 'chansel', 'eye','timedecoding',1, 'nrandperm',100);
  accrand{2,h} = squeeze(mean(accrand{2,h},4));

  
  for n=1:100
    Ir{n}=[];
    for k=1:110
      tmp1 = [];
      tmp2 = [];
      for l=1:5
        tmp1 = [tmp1; datrand{1,h}.testdata_fold{n,l,k}];
        xr{k} = tmp1;
        tmp2 = [tmp2; datrand{2,h}.testdata_fold{n,l,k}];
        yr{k} = tmp2;
      end
      MI_rand{n,h}(k,1) = gcmi_cc(xr{k}, yr{k});
    end
    r_rand(n,h) = corr(MI_rand{n,h}, accrand{1,h}(n,:)');
  end
end
  
save([results_dir, 'MI/', sprintf('sub%02d_mutual_information.mat', subj)], 'accuracy', 'accrand', 'MI', 'MI_rand', 'r', 'r_rand', 'dat', 'datrand')
