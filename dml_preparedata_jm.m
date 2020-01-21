function [traindata, testdata, traindesign, testdesign, P] = dml_preparedata_jm(data, trlnum, tpoint, do_prewhiten, covar)
if ~exist('do_prewhiten', 'var'), do_prewhiten = false; end

ntest  = zeros(numel(data),1);
ntrain = zeros(numel(data),1);
for k=1:numel(data)
  ntrl          = size(data{k}.trial,1);
  testdata{k}   = data{k}.trial(trlnum, :, tpoint);
  traindata{k}  = data{k}.trial(setdiff(1:ntrl, trlnum), :, tpoint);
  
  ntest(k,1)  = size(testdata{k},1);
  ntrain(k,1) = size(traindata{k},1);
end

if do_prewhiten==1 % prewhiten only based on the training set. prewhiten over trials, potentially average over time
  tmp = cat(1,traindata{:});
  M   = mean(tmp);
  tmp = tmp-M(ones(size(tmp,1),1),:);
  covar = tmp'*tmp;
  
  thr     = 0.99;
  [u,s,v] = svd(covar);
  diagS   = diag(s);
  sel     = find(cumsum(diagS)./sum(diagS)<=thr);  
  P       = diag(1./sqrt(diagS(sel)))*u(:,sel)';

  for k=1:numel(traindata)
    traindata{k} = (traindata{k} - M(ones(ntrain(k),1), :))*P';
    testdata{k}  = (testdata{k}  - M(ones(ntest(k), 1), :))*P';
  end
else
  P = [];
end

traindata = cat(1,traindata{:});
testdata  = cat(1,testdata{:});

traindesign = zeros(0,1);
testdesign  = zeros(0, 1);
for k=1:numel(data)
  traindesign = cat(1,traindesign,k.*ones(ntrain(k),1));
  testdesign  = cat(1,testdesign, k.*ones(ntest(k), 1));
end
