function [data, cov, cov_inv] = prewhiten_data(data)
% prewhitens the data by computing the covariance over the entire dataset
% and dividing the data by the average covariance over datasets.
% this function assumes multiple timelockstructures in a cell, with dimord
% 'rpt_chan_time'.

cfg=[];
cfg.removemean = 'no';
cfg.covariance = 'yes';
for k=1:numel(data)
  cov(k) = ft_timelockanalysis(cfg, data{k});
end
cov = mean(cat(3,cov(:).cov),3);
cov_inv = cov^-0.5;

for k=1:numel(data)
  for t = 1:size(data{k}.trial,3)
    data{k}.trial(:,:,t) = data{k}.trial(:,:,t)*cov_inv;
  end
end

