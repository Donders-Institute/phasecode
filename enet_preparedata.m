function [traindata, testdata, traindesign, testdesign] = enet_preparedata(data, trlnum, tpoint)

traindesign = [];
testdesign = [];
traindata = [];
testdata = [];
for k=1:numel(data)
  ntrl = size(data{k}.trial,1);
  traindata = [traindata; reshape(data{k}.trial(setdiff(1:ntrl,trlnum),:,tpoint), ntrl-1, [])];
  testdata = [testdata; reshape(data{k}.trial(trlnum,:,tpoint), numel(trlnum), [])];
  traindesign = [traindesign; k*ones(ntrl-1,1)];
  testdesign = [testdesign; k*ones(numel(trlnum),1)];
end
