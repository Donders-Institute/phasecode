function analysis_dics_ami(subj)

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

cfgf=[];
cfgf.output = 'fourier';
cfgf.method = 'mtmfft';
cfgf.taper = 'hanning';
cfgf.foi = 10;

for k=1:numel(dat)
  data=dat{k};

% if you want a hemisphere specific source, use the following.  
%   cfg=[];
%   cfg.trials = data.trialinfo(:,1)==1; % for left, 2 for right hemifield
% %   attention
%   data = ft_selectdata(cfg, data);
  
  cfg=[];
  cfg.toilim = [-1+1/fs 0];
  cfg.trials = data.trialinfo(:,1)==2;
  dataR = ft_redefinetrial(cfg, data);
  cfg.trials = data.trialinfo(:,1)==1;
  dataL = ft_redefinetrial(cfg, data);
  
  data_all = ft_appenddata(cfg, dataL, dataR);
  
  freq_all = ft_freqanalysis(cfgf, data_all);
  freqL = ft_freqanalysis(cfgf, dataL);
  freqR = ft_freqanalysis(cfgf, dataR);
  
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
  
  sourceL = ft_sourceanalysis(cfg, ft_checkdata(freqL, 'cmbrepresentation', 'fullfast'));
  sourceR   = ft_sourceanalysis(cfg, ft_checkdata(freqR, 'cmbrepresentation', 'fullfast'));
  
  source_inc = sourceL;
  source_inc.avg.pow = (sourceL.avg.pow-sourceR.avg.pow)./(sourceL.avg.pow+sourceR.avg.pow);
  
  
  % projection matrix to get from fourier to power
  nrpt = numel(freqL.cumtapcnt);
  ntap = freqL.cumtapcnt(1);
  
  ix = reshape(repmat(1:nrpt,[ntap 1]),[],1);
  iy = 1:(nrpt*ntap);
  iz = ones(nrpt*ntap,1)./ntap;
  P  = sparse(iy,ix,iz,nrpt*ntap,nrpt);
  
  % Compute single trial power
  F = zeros(numel(source.inside), numel(freqL.label));
  F(source.inside,:) = cat(1, source.avg.filter{:});
  
  nrpt = numel(freqL.cumtapcnt);
  ntap = freqL.cumtapcnt(1);
  ix = reshape(repmat(1:nrpt,[ntap 1]),[],1);
  iy = 1:(nrpt*ntap);
  iz = ones(nrpt*ntap,1)./ntap;
  P  = sparse(iy,ix,iz,nrpt*ntap,nrpt);
  tmppowL{k} = rmfield(source,'avg');
  tmppowL{k}.pow = (abs(F*transpose(freqL.fourierspctrm)).^2)*P;
  
  nrpt = numel(freqR.cumtapcnt);
  ntap = freqR.cumtapcnt(1);
  ix = reshape(repmat(1:nrpt,[ntap 1]),[],1);
  iy = 1:(nrpt*ntap);
  iz = ones(nrpt*ntap,1)./ntap;
  P  = sparse(iy,ix,iz,nrpt*ntap,nrpt);
  tmppowR{k} = rmfield(source, 'avg');
  tmppowR{k}.pow = (abs(F*transpose(freqR.fourierspctrm)).^2)*P;
end

% concatenate sourcepower over sessions
powL = rmfield(tmppowL{1}, 'pow');
powL.pow = [];
powR = rmfield(tmppowR{1}, 'pow');
powR.pow = [];
for k=1:numel(dat)
  powL.pow = [powL.pow tmppowL{k}.pow];
  powR.pow = [powR.pow tmppowR{k}.pow];
end

% Compute T statistics
% nL = size(powL.pow,2);
% nR = size(powR.pow,2);
% cfg = [];
% cfg.ivar = 1;
% design = [ones(1,nL) 2*ones(1,nR);1:nL 1:nR];
% cfg.design = design;
% cfg.method = 'analytic';
% cfg.statistic = 'indepsamplesT';
% Tval = ft_sourcestatistics(cfg, powL, powR);

powAMI = rmfield(powL, 'pow');
powAMI.pow = (mean(powL.pow,2)-mean(powR.pow,2))./(mean(powL.pow,2)+mean(powR.pow,2));

load standard_sourcemodel3d6mm.mat
ixL = find(sourcemodel.pos(:,1)<0);
ixR = find(sourcemodel.pos(:,1)>0);
[Tmax(1), maxidx(1)] = max(Tval.stat(ixL));
[Tmax(2), maxidx(2)] = min(Tval.stat(ixR));
maxidx(1) = ixL(maxidx(1));
maxidx(2) = ixR(maxidx(2));

save([projectdir, sprintf('results/freq/sub%02d_dics_alpha', subj)], 'powAMI', 'maxidx','Tmax', 'Tval', 'powL', 'powR', 'sourcemodel_ses', 'dicsfilter');
