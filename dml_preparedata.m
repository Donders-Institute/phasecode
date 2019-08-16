function [traindata, testdata, traindesign, testdesign] = dml_preparedata(data, trlnum, tpoint, do_prewhiten, cov)
if ~exist('do_prewhiten', 'var'), do_prewhiten = false; end

ntrl = size(data{1}.trial,1);
 
for k=1:numel(data)
  testdata{k} = data{k}.trial(trlnum, :, tpoint);
  traindata{k} = data{k}.trial(setdiff(1:ntrl, trlnum), :, tpoint);
end

if do_prewhiten % prewhiten only based on the training set
  cov = squeeze(mean(cov(setdiff(1:ntrl, trlnum),:,:),1));
  cov_inv = cov^-0.5;
  for k=1:numel(data)
    traindata{k} = traindata{k}*cov_inv;
    testdata{k} = testdata{k}*cov_inv;
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
