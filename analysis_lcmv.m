function analysis_lcmv(subj)

datainfo;

%% load anatomical data
mridir = [projectdir, sprintf('results/anatomy/sub%02d/', subj)];
load([mridir, 'sourcemodel3d'])
load([mridir, 'headmodel'])
headmodel = ft_convert_units(headmodel, 'm');
sourcemodel = ft_convert_units(sourcemodel, 'm');

% index of maximal gamma power increase
load([projectdir, sprintf('results/freq/sub%02d_dics_gamma', subj)], 'maxidx', 'Tval')
if ~isequal(Tval.pos, sourcemodel.pos)
  error('maxidx does not correspond to sourcemodel positions')
end


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
  data = dat{k};
  
  cfg         = [];
  cfg.latency = [-0.25+1/fs 1.2];
  data        = ft_selectdata(cfg, data);
  
  cfg                        = [];
  cfg.preproc.demean         = 'yes';
  cfg.preproc.baselinewindow = [-0.1+1/fs 0];
  cfg.removemean             = 'no';
  cfg.keeptrials             = 'yes';
  tlck                       = ft_timelockanalysis(cfg, data); % use this for virtual channel
  cfg.covariance             = 'yes';
  cfg.keeptrials             = 'no';
  cfg.latency                = [-0.1+1/fs 1];
  tlck_cov                   = ft_timelockanalysis(cfg, data); % use this for lcmv
  
  % prepare anatomical data with session specific leadfields
  cfg                = [];
  cfg.headmodel      = headmodel;
  cfg.grid           = sourcemodel;
  cfg.grad           = ft_convert_units(tlck.grad, 'm');
  cfg.channel        = tlck.label;
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
  
  sesF      = zeros(size(source.pos,1),numel(tlck_cov.label));
  sesF(source.inside,:) = cat(1,source.avg.filter{:});
  F{k} = sesF;
  
  sesL      = zeros(numel(tlck_cov.label), size(source.pos,1));
  sesL(:,source.inside) = cat(2,source.leadfield{:});
  L{k} = sesL;
  
  tmpvirtualchan = keepfields(tlck, {'time', 'dimord','trialinfo', 'sampleinfo'});
  tmpvirtualchan.label{1} = 'maxgamma';
  tmpvirtualchan.trial = zeros(size(tlck.trial,1),1,size(tlck.trial,3));
  for t = 1:size(tlck.trial,3)
    tmpvirtualchan.trial(:,1,t) = transpose(sesF(maxidx, :)*tlck.trial(:,:,t)');
  end
  
  if k==1
    virtualchan = tmpvirtualchan;
  else
    virtualchan.trial = [virtualchan.trial; tmpvirtualchan.trial];
    virtualchan.trialinfo = [virtualchan.trialinfo; tmpvirtualchan.trialinfo];
    virtualchan.sampleinfo = [virtualchan.sampleinfo; tmpvirtualchan.sampleinfo];    
  end 
end

save([projectdir, sprintf('results/tlck/sub%02d_virtualchan', subj)], 'virtualchan', 'F', 'L')







