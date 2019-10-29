function analysis_dics(subj)

% use a DICS beamformer to find the location with the highest gamma power
% increase from baseline. Do this seperately for each session, as they
% have different grad structures.

datainfo;

%% load anatomical data
mridir = [projectdir, sprintf('results/anatomy/sub%02d/', subj)];
load([mridir, 'sourcemodel3d'])
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
load([projectdir, sprintf('results/freq/sub%02d_peakfreq.mat', subj)], 'gamma');
f = gamma.peakfreq;

cfgf=[];
cfgf.output = 'fourier';
cfgf.method = 'mtmfft';
cfgf.taper = 'hanning';
cfgf.foi = f;

for k=1:numel(dat)
  data=dat{k};

% if you want a hemisphere specific source, use the following.  
%   cfg=[];
%   cfg.trials = data.trialinfo(:,1)==1; % for left, 2 for right hemifield
% %   attention
%   data = ft_selectdata(cfg, data);
  
  cfg=[];
  cfg.toilim = [-1+1/fs 0];
  data_bl = ft_redefinetrial(cfg, data);
  cfg.toilim = [1/fs 1];
  data_stim = ft_redefinetrial(cfg, data);
  
  data_all = ft_appenddata(cfg, data_stim, data_bl);
  
  freq_all = ft_freqanalysis(cfgf, data_all);
  freq_stim = ft_freqanalysis(cfgf, data_stim);
  freq_bl = ft_freqanalysis(cfgf, data_bl);
  
  % prepare anatomical data with session specific leadfields
  cfg               = [];
  cfg.headmodel     = headmodel;
  cfg.grid          = sourcemodel;
  cfg.grad          = ft_convert_units(freq_all.grad, 'm');
  cfg.channel       = freq_all.label;
  sourcemodel_ses{k} = ft_prepare_leadfield(cfg);
  
  cfg                 = [];
  cfg.method          = 'dics';
  cfg.frequency       = f;
  cfg.headmodel       = headmodel;
  cfg.sourcemodel     = sourcemodel_ses{k};
  cfg.dics.keepfilter = 'yes';
  cfg.dics.fixedori   = 'yes';
  cfg.dics.realfilter = 'yes';
  cfg.dics.lambda     = '100%';
  source              = ft_sourceanalysis(cfg, ft_checkdata(freq_all, 'cmbrepresentation', 'fullfast'));
  dicsfilter{k} = source.avg.filter;
  cfg.sourcemodel.filter = dicsfilter{k};
  
  source_stim = ft_sourceanalysis(cfg, ft_checkdata(freq_stim, 'cmbrepresentation', 'fullfast'));
  source_bl   = ft_sourceanalysis(cfg, ft_checkdata(freq_bl, 'cmbrepresentation', 'fullfast'));
  
  source_inc = source_stim;
  source_inc.avg.pow = source_stim.avg.pow./source_bl.avg.pow-1;
  
  
  % projection matrix to get from fourier to power
  nrpt = numel(freq_stim.cumtapcnt);
  ntap = freq_stim.cumtapcnt(1);
  
  ix = reshape(repmat(1:nrpt,[ntap 1]),[],1);
  iy = 1:(nrpt*ntap);
  iz = ones(nrpt*ntap,1)./ntap;
  P  = sparse(iy,ix,iz,nrpt*ntap,nrpt);
  
  % Compute single trial power
  F = zeros(numel(source.inside), numel(freq_stim.label));
  F(source.inside,:) = cat(1, source.avg.filter{:});
  
  tmppow_stim{k} = rmfield(source,'avg');
  tmppow_stim{k}.pow = (abs(F*transpose(freq_stim.fourierspctrm)).^2)*P;
  tmppow_bl{k} = rmfield(source, 'avg');
  tmppow_bl{k}.pow = (abs(F*transpose(freq_bl.fourierspctrm)).^2)*P;
end

% concatenate sourcepower over sessions
pow_stim = rmfield(tmppow_stim{1}, 'pow');
pow_stim.pow = [];
pow_bl = rmfield(tmppow_bl{1}, 'pow');
pow_bl.pow = [];
for k=1:numel(dat)
  pow_stim.pow = [pow_stim.pow tmppow_stim{k}.pow];
  pow_bl.pow = [pow_bl.pow tmppow_bl{k}.pow];
end

% Compute T statistics
n = size(pow_stim.pow,2);
cfg = [];
cfg.ivar = 1;
cfg.uvar = 2;
design = [ones(1,n) 2*ones(1,n);1:n 1:n];
cfg.design = design;
cfg.method = 'analytic';
cfg.statistic = 'depsamplesT';
Tval = ft_sourcestatistics(cfg, pow_stim, pow_bl);

[~, maxidx] = max(Tval.stat);

save([projectdir, sprintf('results/freq/sub%02d_dics_gamma', subj)], 'maxidx', 'Tval', 'pow_stim', 'pow_bl', 'sourcemodel_ses', 'dicsfilter');
