 function [stat, cfg] = statfun_perm(cfg, dat, randdat)

% STATFUN_xxx is a function for computing a statistic for the relation
% between biological data and a design vector containing trial
% classifications or another independent variable
%
% This function is called by STATISTICS_MONTECARLO, where you can specify
% cfg.statistic = 'xxx' which will be evaluated as statfun_xxx.
%
% The external interface of this function has to be
%   [s] = statfun_xxx(cfg, dat, design);
% where
%   dat    contains the biological data, Nvoxels x Nreplications
%   design contains the independent variable,  1 x Nreplications
%
% Additional settings can be passed through to this function using
% the cfg structure.
%
% STATFUN_COSINEFIT fits a cosine to the data dat. the independent
% variable design should contain angular values, bounded by -pi and pi.
%
% the output s is a structure containing the statistic, as specified by
%  cfg.cosinefit.statistic. this can either be the amplitude (default) of the fit, angle   (giving the preferred angle), complex (giving both angle and amplitude in a complex number), or fit (giving the percentage of explained variance).
% additional fields in the output=structure are:
%  s.r      percentage of explained variance.
%  s.offset the DC-component of the fit.
%
% Additional cfg-options are:
%  cfg.cosinefit.phi = [] (default) or angular value between -pi and pi, estimates the amplitude of a cosine at a given angle.

% Copyright (C) 2006 Jan-Mathijs Schoffelen

if ~isfield(cfg, 'numrandomization'), numrandomization = 1; end
if ~isfield(cfg, 'ivar'), error('ivar was not provided'); end
if ~isfield(cfg, 'uvar'), error('uvar was not provided'); end
cfg.alpha            = ft_getopt(cfg, 'alpha',      0.05);
cfg.tail             = ft_getopt(cfg, 'tail',       0);
cfg.correctm         = ft_getopt(cfg, 'correctm',   'no');
cfg.resampling       = ft_getopt(cfg, 'resampling', 'permutation');
cfg.clusteralpha     = ft_getopt(cfg, 'clusteralpha',      0.05);
cfg.clustertail      = ft_getopt(cfg, 'clustertail',      0);
cfg.clusterstatistic = ft_getopt(cfg, 'clusterstatistic', 'maxsum');
cfg.clusterthreshold = ft_getopt(cfg, 'clusterthreshold', 'parametric');
cfg.correctm         = ft_getopt(cfg, 'correctm', []);

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

uncorrected_p = sum(statobs>statrand,2)./Nrand;
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

