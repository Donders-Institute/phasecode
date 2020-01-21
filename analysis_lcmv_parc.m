function [source_parc] = analysis_lcmv_parc(subj)
% Model time courses on the source level (cortical sheet) using an LCMV 
% beamformer. Dipoles are combined into anatomically defined
% parcels.
%
datainfo;

%% load data
mridir = [projectdir, sprintf('results/anatomy/sub%02d/', subj)];
load([mridir, 'sourcemodel2d'])
load([mridir, 'headmodel'])
headmodel = ft_convert_units(headmodel, 'm');
sourcemodel = ft_convert_units(sourcemodel, 'm');


cnt=1;
for ses=subjects(subj).validsessions
  filename = [datadir, sprintf('sub%02d/meg%02d/sub%02d-meg%02d_cleandata.mat', subj, ses, subj, ses)];
  dat{cnt} = load(filename, 'data');
  dat{cnt} = dat{cnt}.data;
  cnt=cnt+1;
end

for k=1:numel(dat)
  data = dat{k};
  
  cfg         = [];
  cfg.latency = [-0.25+1/fs 1.2];
  data        = ft_selectdata(cfg, data);
  cfg                        = [];
  cfg.preproc.demean         = 'yes';
  cfg.preproc.baselinewindow = [-0.1+1/fs 0];
  cfg.removemean             = 'no';
  cfg.keeptrials             = 'yes';
  tlck{k}                       = ft_timelockanalysis(cfg, data); % use this for virtual channel
  
  cfg.covariance             = 'yes';
  cfg.keeptrials             = 'no';
  cfg.latency                = [-0.1+1/fs 1];
  tlck_cov                   = ft_timelockanalysis(cfg, data); % use this for lcmv
  
  cfg                = [];
  cfg.headmodel      = headmodel;
  cfg.grid           = sourcemodel;
  cfg.grad           = ft_convert_units(tlck{k}.grad, 'm');
  cfg.channel        = tlck{k}.label;
  sourcemodel_ses{k} = ft_prepare_leadfield(cfg);
  
  
  
cfg = [];
cfg.method = 'lcmv';
cfg.headmodel = headmodel;
cfg.grid      = sourcemodel;
cfg.keepleadfield = 'yes';
cfg.lcmv.keepfilter = 'yes';
cfg.lcmv.fixedori   = 'yes';
cfg.lcmv.lambda     = '100%';
cfg.lcmv.keepleadfield = 'yes';
cfg.lcmv.keepori = 'yes';
cfg.lcmv.weightnorm = 'unitnoisegain';
source = ft_sourceanalysis(cfg, tlck);

F      = zeros(size(source.pos,1),numel(tlck.label));
F(source.inside,:) = cat(1,source.avg.filter{:});

L      = zeros(numel(tlck.label), size(source.pos,1));
L(:,source.inside) = cat(2,source.leadfield{:});

  sesF      = zeros(size(source.pos,1), numel(tlck_cov.label));
  sesF(source.inside,:) = cat(1,source.avg.filter{:});
  F{k} = sesF;
  
  sesL      = zeros(numel(tlck_cov.label), size(source.pos,1));
  sesL(:,source.inside) = cat(2,source.leadfield{:});
  L{k} = sesL;
end



% prepare the cfg for pca
cfg                       = [];
cfg.method                = 'pca';

for k = 1:numel(selparc)
    tmpF = F(atlas.parcellation==selparc(k),:);
    tmp.trial = tmpF*data.trial;
    for l=1:size(tmpF,1)
        tmplabel{l} = sprintf('loc%03d', l);
    end
    tmp.label = tmplabel;
    clear tmplabel
    cfg.comment = sprintf('Select the spatial filters for all nodes belonging to parcel %s. Create the source time courses for the nodes in this parcel by multiplying the spatial filters with the channel level data. Give this as input to ft_componentanalysis.', atlas.parcellationlabel{selparc(k)});
    tmpcomp   = ft_componentanalysis(cfg, tmp);

    tmpL = L(:,atlas.parcellation==selparc(k));
    
    source_parc.F{k}     = tmpcomp.unmixing*tmpF;
    source_parc.L{k}     = tmpL*tmpcomp.unmixing';
    source_parc.avg(k,:) = source_parc.F{k}(1,:)*tlck.avg;
    
    source_parc.F{k} = source_parc.F{k}(1,:);
end


source_parc.F = cat(1, source_parc.F{:});



