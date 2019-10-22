function [data, cov, cov_inv] = prewhiten_data(data)
% prewhitens the data by computing the covariance over the entire dataset
% and dividing the data by the average covariance over datasets.
% this function assumes multiple timelockstructures in a cell, with dimord
% 'rpt_chan_time'.

cfg=[];
cfg.removemean = 'no';
cfg.covariance = 'yes';
cfg.keeptrials = 'yes';
for k=1:numel(data)
  cov{k} = timelockprewhiten(data{k}); % this line replaces the following to increase computation speed. 
  %{
  cov{k} = ft_timelockanalysis(cfg, data{k});
  cov{k} = cov{k}.cov;
  %}
end
cov = mean(cat(4,cov{:}),4);
cov_avg= squeeze(mean(cov,1)); % FIXME implement a shrinkage transform (Ledoit and Wolf, 2004), to prevent rank-deficiency
cov_inv = cov_avg^-0.5;

for k=1:numel(data)
  for t = 1:size(data{k}.trial,3)
    data{k}.trial(:,:,t) = data{k}.trial(:,:,t)*cov_inv;
  end
end

