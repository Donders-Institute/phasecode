function [traindata, testdata, traindesign, testdesign, covar] = dml_preparedata(data, trlnum, tpoint, do_prewhiten, covar)
if ~exist('do_prewhiten', 'var'), do_prewhiten = false; end

ntrl = size(data{1}.trial,1);
 
for k=1:numel(data)
  testdata{k} = data{k}.trial(trlnum, :, tpoint);
  data{k}.trial = data{k}.trial(setdiff(1:ntrl, trlnum), :, :);
  traindata{k} = data{k}.trial(:,:,tpoint);
end

if do_prewhiten % prewhiten only based on the training set. prewhiten over trials, potentially average over time
  if ~exist('cov', 'var') || isempty(covar)
    for k=1:numel(data)
      data{k}.trial = permute(data{k}.trial(:,:,:), [3 2 1]);
      data{k}.time = 1:size(data{k}.trial,3);
    end
    [~, covar] = prewhiten_data(data);
    covar = squeeze(mean(covar,1));
  end
  % instead of inversing the covariance matrix (which could run into
  % numerical problems if it's rank deficient), first shrink the number of
  % components (and thus svm features).
%   covar_inv = covar^-0.5;
  [u,s,v]=svd(covar);
  diagS=diag(s);
  sel=find(cumsum(diagS)./sum(diagS)<=0.99);  
  P=diag(1./sqrt(diagS(sel)))*u(:,sel)';

  for k=1:numel(data)
        traindata{k} = traindata{k}*P';
    testdata{k} = testdata{k}*P';
%     traindata{k} = traindata{k}*covar_inv;
%     testdata{k} = testdata{k}*covar_inv;
  end
end

traindata = cat(1,traindata{:});
testdata = cat(1,testdata{:});
traindesign = [];
testdesign = [];
for k=1:numel(data)
  traindesign = [traindesign; k*ones(ntrl-numel(trlnum),1)];
  testdesign = [testdesign; k*ones(numel(trlnum),1)];
end
