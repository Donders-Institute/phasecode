function [stat, cfg, dat] = statfun_yuenTtest(cfg, dat, design)

cfg = struct(cfg);

if ~isfield(cfg, 'yuen'),         cfg.yuen         = [];  end
if ~isfield(cfg.yuen, 'percent'), cfg.yuen.percent = 0.2; end %default handling is done by low level function
if ~isfield(cfg, 'preconditionflag'), cfg.preconditionflag = 0; end

cfg.yuen.type = ft_getopt(cfg.yuen, 'type', 'indepsamples'); 
if isempty(exist('yuent', 'file')),
  addpath('/home/language/jansch/toolboxes/robuststats');
end
if ~isfield(cfg, 'ivar'), cfg.ivar = 1; end

if cfg.preconditionflag
  fprintf('removing the marginal means from the individual cells, updated data will be used for further computations\n');
  cell_id = unique(design(cfg.ivar,:));
  for k = 1:length(cell_id)
    tmpsel = design(cfg.ivar,:)==cell_id(k);
    dat(:,tmpsel) = dat(:,tmpsel) - repmat(mean(dat(:,tmpsel),2),[1 sum(tmpsel)]);
  end
  stat = [];
end

sel1 = find(design(cfg.ivar,:)==1);
sel2 = find(design(cfg.ivar,:)==2);

switch cfg.yuen.type
  case 'indepsamples'
    stat.stat = yuent(dat(:,sel1), dat(:,sel2), cfg.yuen.percent, 2);
  case 'depsamples'
    uvar1 = design(cfg.uvar,sel1);
    uvar2 = design(cfg.uvar,sel2);
    
    [~, i1,i2] = intersect(uvar1, uvar2);
    
    stat.stat = yuent_dep(dat(:,sel1(i1)), dat(:,sel2(i2)), cfg.yuen.percent, 2);
  otherwise
end