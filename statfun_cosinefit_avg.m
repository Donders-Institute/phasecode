function s = statfun_cosinefit_avg(cfg, dat, design)

if ~isfield(cfg, 'cosinefit'),             cfg.cosinefit = [];                      end
if ~isfield(cfg.cosinefit, 'statistic'),   cfg.cosinefit.statistic   = 'amplitude'; end
if ~isfield(cfg.cosinefit, 'avgoverchan'), cfg.cosinefit.avgoverchan = 'no';        end
if ~isfield(cfg, 'dim'),                   cfg.dim = size(design);                  end

tmpcfg.cosinefit           = cfg.cosinefit;
tmpcfg.cosinefit.statistic = 'complex'; %always take complex output

if size(dat,3)~=1,                  error('size(dat,3) should be 1'); end
if size(dat,1)~=prod(size(design(cfg.factor,:))), error('dimensionality of the design does not correspond with the dimensionality of the data'); end

dim     = cfg.dim;
if length(cfg.dim) > 2,
  design = reshape(design, [size(design,1) dim(1)*dim(2)*dim(3)]);
end

tmp       = reshape(nanmean(dat,2), [dim(1) dim(2)*dim(3)]);
tmpdesign = reshape(design(cfg.factor,:), [dim(1) dim(2)*dim(3)]);
s         = statfun_cosinefit(tmpcfg, tmp', tmpdesign');
s.r     = reshape(s.r   ,   [dim(2) dim(3)]);
s.s     = reshape(s.s   ,   [dim(2) dim(3)]);
s.offset= reshape(s.offset, [dim(2) dim(3)]);
s.phs   = design;
s.dat   = reshape(tmp,    [dim(1) dim(2) dim(3)]);
s.stat  = reshape(s.stat, [dim(2) dim(3)]);

if strcmp(cfg.cosinefit.statistic, 'tstat'),
  s.stat = s.stat./s.s;
end

if strcmp(cfg.cosinefit.avgoverchan, 'complex'),
  s.stat = mean(s.stat);
elseif strcmp(cfg.cosinefit.avgoverchan, 'abs'),
  s.stat = mean(abs(s.stat));
end

if strcmp(cfg.cosinefit.statistic, 'amplitude'),
  s.phi   = angle(s.stat);
  s.stat  = abs(s.stat);
elseif strcmp(cfg.cosinefit.statistic, 'complex'),
  s.stat  = s.stat;
elseif strcmp(cfg.cosinefit.statistic, 'tstat'),
  s.phi   = angle(s.stat);
  s.stat  = abs(s.stat);
else
  error('cfg.cosinefit.statistic should either be amplitude or complex');
end

s.cfg = cfg;
