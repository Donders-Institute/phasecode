 function [stat, cfg] = statfun_perm(cfg, dat, randdat)

% STATFUN_PERM is a function for computing a statistic based on a variable
% and a precpmputed permutation distribution. Only working for cfg.tail =1;
%
% Copyright (C) 2020 Mats van Es
warning('this function currently only works with cfg.tail=1')
if ~isfield(cfg, 'numrandomization'), numrandomization = 1; end
if ~isfield(cfg, 'ivar'), error('ivar was not provided'); end
if ~isfield(cfg, 'uvar'), error('uvar was not provided'); end
cfg.alpha            = ft_getopt(cfg, 'alpha',      0.05);
cfg.tail             = ft_getopt(cfg, 'tail',       1);
cfg.correctm         = ft_getopt(cfg, 'correctm',   'no');
cfg.resampling       = ft_getopt(cfg, 'resampling', 'permutation');
cfg.clusteralpha     = ft_getopt(cfg, 'clusteralpha',      0.05);
cfg.clustertail      = ft_getopt(cfg, 'clustertail',      0);
cfg.clusterstatistic = ft_getopt(cfg, 'clusterstatistic', 'maxsum');
cfg.clusterthreshold = ft_getopt(cfg, 'clusterthreshold', 'parametric');
cfg.correctm         = ft_getopt(cfg, 'correctm', []);
if cfg.tail~=1
    error('this function only works with cfg.tail=1');
end
Nrand                = getfield(cfg, 'numrandomization');
uvar                 = getfield(cfg, 'uvar');
ivar                 = getfield(cfg, 'ivar');


%---create and check independent variable
if size(randdat, 2) ~= size(dat,2), error('random data is incompatible with input-data'); end
% if size(randdat, 1) ~= size(dat,1) && size(randdat,1)==1,
%   quickflag = 1;
%   repdim    = [size(dat,1) 1];
% else
%   quickflag = 0;
%   repdim    = [1 1];
% end  

statobs = nanmean(dat,uvar); % should this be normalized?

Nperm = size(randdat, setdiff(1:ndims(randdat), [uvar, ivar]));
Nrepl = size(dat, uvar);
for i=1:Nrand
  for j=1:Nrepl
      resample(i,j) = randi(Nperm,1);
  end
end
    
for i=1:Nrand
  for j=1:Nrepl
    tmp(:,j) = randdat(:,j, resample(i,j));
  end
  statrand(:,i) = nanmean(tmp,2);
end

uncorrected_p = sum(statobs<statrand,2)./Nrand;
uncorrected_p = max(uncorrected_p, 1./Nrand); % minimum p-value depends on number of randomizations
if ~isfield(cfg, 'dim')
  cfg.dim = size(statobs);
end
if strcmp(cfg.correctm, 'cluster');
  cfg.feedback = 'no';
  
  stat = clusterstat(cfg, statrand, statobs);
  
  stat.stat = reshape(statobs, cfg.dim);
  stat.prob = reshape(stat.prob, cfg.dim);
  stat.uncorrected_p = reshape(uncorrected_p, cfg.dim);
  try
    if cfg.tail>=0
      stat.posclusterslabelmat = reshape(stat.posclusterslabelmat, cfg.dim);
    end
    if cfg.tail<=0
      stat.negclusterslabelmat = reshape(stat.negclusterslabelmat, cfg.dim);
    end
  end
else
  stat = [];
  stat.stat = statobs;
  stat.uncorrectec_p = uncorrected_p;
end

