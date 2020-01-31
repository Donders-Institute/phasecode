function [dat, tlck, F, L] = analysis_lcmv(subj, model)
if nargin<2 || ~exist('model', 'var'), model = '3d'; end
datainfo;

%% load anatomical data
mridir = [projectdir, sprintf('results/anatomy/sub%02d/', subj)];
if strcmp(model, '2d')
  load([mridir, 'sourcemodel2d'])
  load atlas_subparc374_8k.mat
else
  load([mridir, 'sourcemodel3d'])
end
load([mridir, 'headmodel'])
headmodel = ft_convert_units(headmodel, 'm');
sourcemodel = ft_convert_units(sourcemodel, 'm');


%% load MEG data
cnt=1;
for ses=subjects(subj).validsessions
  filename = [datadir, sprintf('sub%02d/meg%02d/sub%02d-meg%02d_cleandata.mat', subj, ses, subj, ses)];
  dat{cnt} = load(filename, 'data');
  dat{cnt} = dat{cnt}.data;
  cnt=cnt+1;
end

% prepare configurations and variables
fs = dat{1}.fsample;


for k=1:numel(dat)
  cfg         = [];
  cfg.latency = [-0.25+1/fs 1.2];
  dat{k}        = ft_selectdata(cfg, dat{k});
  
  cfg                        = [];
  cfg.preproc.demean         = 'yes';
  cfg.preproc.baselinewindow = [-0.1+1/fs 0];
  cfg.removemean             = 'no';
  cfg.keeptrials             = 'yes';
  tlck{k}                       = ft_timelockanalysis(cfg, dat{k}); % use this for virtual channel
  cfg.covariance             = 'yes';
  cfg.keeptrials             = 'no';
  cfg.latency                = [-0.1+1/fs 1];
  tlck_cov                   = ft_timelockanalysis(cfg, dat{k}); % use this for lcmv
  
  % prepare anatomical data with session specific leadfields
  cfg                = [];
  cfg.headmodel      = headmodel;
  cfg.grid           = sourcemodel;
  cfg.grad           = ft_convert_units(tlck{k}.grad, 'm');
  cfg.channel        = tlck{k}.label;
  sourcemodel_ses{k} = ft_prepare_leadfield(cfg);
  
  cfg                    = [];
  cfg.method             = 'lcmv';
  cfg.headmodel          = headmodel;
  cfg.sourcemodel        = sourcemodel_ses{k};
  cfg.keepleadfield      = 'yes';
  cfg.lcmv.keepfilter    = 'yes';
  cfg.lcmv.fixedori      = 'yes';
  cfg.lcmv.lambda        = '100%';
  cfg.lcmv.keepleadfield = 'yes';
  cfg.lcmv.keepori       = 'yes';
  cfg.lcmv.weightnorm    = 'unitnoisegain';
  source                 = ft_sourceanalysis(cfg, tlck_cov);
  
  sesF      = zeros(size(source.pos,1), numel(tlck_cov.label));
  sesF(source.inside,:) = cat(1,source.avg.filter{:});
  F{k} = sesF;
  
  sesL      = zeros(numel(tlck_cov.label), size(source.pos,1));
  sesL(:,source.inside) = cat(2,source.leadfield{:});
  L{k} = sesL;
end



