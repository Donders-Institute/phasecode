function [phasebin, phase, dist, time] = analysis_alphaphase(data, bpfreq, centerphase, latency)
% This function load the original data and for each trial finds the
% occipital channel that represents the alpha phase time course of the
% occipital average best. For every time point, it then divides the trials
% into four phase bins.

data_orig=data;

cfg=[];
cfg.channel = 'MEG';%{'MLO', 'MZO', 'MRO'};
cfg.latency = latency;
data = ft_selectdata(cfg, data_orig);

% bandpass filter in alpha band and hilbert transform
cfg=[];
cfg.bpfilter = 'yes';
cfg.bpfreq = bpfreq;
cfg.bpfilttype = 'firws';
cfg.hilbert = 'complex';
cfg.usefftfilt = 'yes';
data_alpha = ft_preprocessing(cfg, data);
alphacomplex = cat(3,data_alpha.trial{:});
alphaprincipal = zeros(size(alphacomplex,3), size(alphacomplex,2));

% define the phase based on an SVD on the complex value over channels
for k=1:size(alphacomplex,3)
    [u,s,v]=svd(alphacomplex(:,:,k),'econ');
    alphaprincipal(k,:) = u(:,1)'*alphacomplex(:,:,k);% make sure phases are between 0 and 2*pi
end
phase = angle(alphaprincipal)+pi;

% for each time window, select only those trials in which alpha phase is in
% the same phase bin.

% define time axis
time = data_alpha.time{1};
fs = data_alpha.fsample;
dist = zeros(size(phase));
twindow = 0.02;
nsample = twindow*fs;

tmpcenterphase = [centerphase 2*pi];
for k=1:numel(time)
    ptmp = phase(:,k);
    [tmpdist, idx] = min(transpose(abs(ptmp-tmpcenterphase))); % distance to center phases.
    idx(idx==max(idx)) = min(idx);
    dist(:,k) = tmpdist;
    phasebin(:,k) = idx;
end


