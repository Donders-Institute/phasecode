function [cov] = timelockprewhiten(data)
% This function copies the essential steps from ft_timelockanalysis,
% without all the administrative code. This function is aimed to speed up
% the processing pipeline of phasecode's analysis_decoding.m (through
% dml_preparedata.m and prewhiten_data.m).
datacov=data;
[nrpt, nchan, nsmp] = size(datacov.trial);
cov = nan(nrpt, nchan, nchan);

for k=1:nrpt
  dat    = reshape(datacov.trial(k,:,:), [nchan nsmp]);
  datsmp = isfinite(dat);
  numsmp = sum(datsmp(1,:));
  dat(~datsmp)  = 0;
  cov(k,:,:) = dat*dat'./numsmp;
end