function [traindata, testdata, traindesign, testdesign, cov] = dml_preparedata(data, trlnum, tpoint, do_prewhiten, cov)
if ~exist('do_prewhiten', 'var'), do_prewhiten = false; end

ntrl = size(data{1}.trial,1);
 
for k=1:numel(data)
  testdata{k} = data{k}.trial(trlnum, :, tpoint);
  data{k}.trial = data{k}.trial(setdiff(1:ntrl, trlnum), :, :);
  traindata{k} = data{k}.trial(:,:,tpoint);
end

if do_prewhiten % prewhiten only based on the training set. prewhiten over trials, potentially average over time
  if ~exist('cov', 'var') || isempty(cov)
    for k=1:numel(data)
      data{k}.trial = permute(data{k}.trial(:,:,:), [3 2 1]);
      data{k}.time = 1:size(data{k}.trial,3);
    end
    [~, cov] = prewhiten_data(data);
    cov = squeeze(mean(cov,1));
  end
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
