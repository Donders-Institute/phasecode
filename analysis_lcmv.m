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
x = subjects(subj).sessions;
for k=1:numel(x)
  ses = subjects(subj).sessions(k);
  filename = [datadir, sprintf('3016045.07_matves_%03d_%03d', subj, ses), '/cleandata.mat'];
  filename = [projectdir, sprintf('data/sub%02d/meg%02d/sub%02d-meg%02d_cleandata.mat', subj, ses, subj, ses)];
  dat{k} = load(filename, 'data');
  dat{k} = dat{k}.data;
end

% concatenate those datasets that belong to the same sessions
if numel(x)>3 % MEG system had to be rebooted within 1 session
  for k=1:numel(x)
    tmp = num2str(x(k));
    x(k) = str2num(tmp(1));
  end
  for k=1:3
    idx = find(x==k);
    tmpdat{k} = ft_appenddata([], dat{idx});
    tmpdat{k}.grad = dat{idx(1)}.grad;
  end
  dat = tmpdat;
end


% prepare configurations and variables
fs = dat{1}.fsample;


for k=1:3
  data = dat{k};
  
  cfg         = [];
  cfg.latency = [-0.25+1/fs 1];
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
  cfg               = [];
  cfg.headmodel     = headmodel;
  cfg.grid          = sourcemodel;
  cfg.grad          = ft_convert_units(tlck.grad, 'm');
  cfg.channel       = tlck.label;
  tmpsourcemodel{k} = ft_prepare_leadfield(cfg);
  
  cfg = [];
  cfg.method = 'lcmv';
  cfg.headmodel = headmodel;
  cfg.grid      = tmpsourcemodel{k};
  cfg.keepleadfield = 'yes';
  cfg.lcmv.keepfilter = 'yes';
  cfg.lcmv.fixedori   = 'yes';
  cfg.lcmv.lambda     = '100%';
  cfg.lcmv.keepleadfield = 'yes';
  cfg.lcmv.keepori = 'yes';
  cfg.lcmv.weightnorm = 'unitnoisegain';
  source = ft_sourceanalysis(cfg, tlck_cov);
  
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
  end 
end

save([projectdir, sprintf('results/tlck/sub%02d_virtualchan', subj)], 'virtualchan', 'F', 'L')







