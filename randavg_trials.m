function [data] = randavg_trials(data, ngroups, groupsize)
% this function averages trials without replacements into ngroups of
% 'groupsize' trials. Trials will be randomized according to randum
% beforehand. This function assumes a timelock structure with dimord
% 'rpt_chan_time'.
[~,s2,s3] = size(data.trial);
tmpdata = zeros(ngroups, s2, s3);

idx=1;
for k=1:ngroups
  tmpdata(k,:,:) = mean(data.trial(idx:idx+groupsize-1,:,:),1);
  idx=idx+groupsize;
end

data.trial = tmpdata;

