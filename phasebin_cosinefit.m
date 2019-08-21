
freqs = [4:1:30 32:2:80];
numrand = 50;
subj=4;
centerphase = [0 1/3 2/3 1 4/3 5/3]*pi-pi;


cnt=1;
for f=freqs
  for ii = 1:numrand
    tmp{cnt,ii} = load(sprintf('sub%02d_decoding_f%d_%d', subj, f, ii));
    tmp{cnt,ii} = mean(tmp{cnt,ii}.accuracy,2);
  end
  cnt=cnt+1;
end

for k=1:size(tmp,1)
  for ii=1:numrand
    A(ii,:,k) = tmp{k,ii};
  end
end

dat=[];
dat.dimord = 'rpt_chan_freq';
dat.label{1} = 'accuracy';
dat.freq = freqs;
dat.trial(:,1,:) = reshape(A, [],1, numel(freqs));

nc = numel(centerphase);
design=[];
for k=1:nc
  design = [design; centerphase(k)*ones(numrand,1)];
end


cfg = [];
cfg.method = 'montecarlo';
cfg.statistic = 'cosinefit';
cfg.design = design;
cfg.numrandomization = 100; 
cfg.correctm = 'cluster';
cfg.parameter = 'trial';
stat = ft_freqstatistics(cfg, dat);
% for f=1:nfreq
%   dat.trial(:,1) = reshape(B(:,:,f), [],1);
%   stat(f) = ft_timelockstatistics(cfg, dat);
% end
