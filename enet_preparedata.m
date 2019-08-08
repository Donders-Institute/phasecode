function [traindata, testdata, traindesign, testdesign] = enet_preparedata(data, trlnum, tpoint, do_prewhiten, cov)
if ~exist('do_prewhiten', 'var'), do_prewhiten = false; end

ntrl = size(data{1}.trial,1);
 
cfg=[];
for k=1:numel(data)
  cfg.trials = trlnum;
  cfg.latency = [data{k}.time(tpoint) data{k}.time(tpoint)];
  testdata{k} = ft_selectdata(cfg, data{k});
  testdata{k} = testdata{k}.trial;
  
  cfg.trials = setdiff(1:ntrl, trlnum);
  traindata{k} = ft_selectdata(cfg, data{k});
  traindata{k} = traindata{k}.trial;
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
