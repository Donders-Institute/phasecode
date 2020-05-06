function analysis_timfreq(subj, varargin)

dosave = ft_getopt(varargin, 'dosave', true);

datainfo;
cfg=[];
cnt=1;
for ses=subjects(subj).validsessions
  filename = [datadir, sprintf('sub%02d/meg%02d/sub%02d-meg%02d_cleandata.mat', subj, ses, subj, ses)];
  dat{cnt} = load(filename, 'data');
  dat{cnt} = removefields(dat{cnt}.data, 'elec');
  cnt=cnt+1;
end


% select correct trials
for k=1:numel(dat)
  cfg=[];
  cfg.trials = (dat{k}.trialinfo(:,7)==1) & (dat{k}.trialinfo(:,1)==dat{k}.trialinfo(:,4));
  dat{k} = ft_selectdata(cfg, dat{k});
end


% TFR for each session seperately
for k=1:numel(dat)
  fs = dat{k}.fsample;
  
  % transform to planar gradients
  cfg                 = [];
  cfg.method          = 'template';
  cfg.template        = 'CTF275_neighb.mat';
  cfg.neighbours      = ft_prepare_neighbours(cfg, dat{k});
  cfg.method          = 'sincos';
  dat_pl              = ft_megplanar(cfg, dat{k});
  

  % TFR
  % low frequencies (<30 Hz)
  cfg              = [];
  cfg.output       = 'pow';
  cfg.method       = 'mtmconvol';
  cfg.taper        = 'hanning';
  cfg.foi          = 2:2:30;                         % analysis 2 to 30 Hz in steps of 2 Hz
  cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
  cfg.toi          = -1+1./fs:0.1:1;                  % time window "slides" from -0.8 to 1.25 sec in steps of 0.05 sec (50 ms)
  cfg.keeptrials   = 'yes';
  cfg.pad          = 4;
  tfr_low         = ft_freqanalysis(cfg, dat_pl);
  
  cfg=[];
  tfr_plcmb_low{k} = ft_combineplanar(cfg, tfr_low);
  
  % high frequencies
  cfg = [];
  cfg.output       = 'pow';
  cfg.method       = 'mtmconvol';
  cfg.foi          = 30:4:80;
  cfg.t_ftimwin    = 5./cfg.foi;
  cfg.tapsmofrq    = 0.2 *cfg.foi; %width of frequency smoothing
  cfg.toi          = -1+1./fs:0.05:1;
  cfg.keeptrials   = 'yes';
  cfg.pad          = 4;
  tfr_high         = ft_freqanalysis(cfg, dat_pl);
  cfg=[];
  tfr_plcmb_high{k} = ft_combineplanar([], tfr_high);
end

cfg=[];
cfg.appenddim = 'rpt';
cfg.parameter = 'powspctrm';
tfr_plcmb_low = ft_appendfreq(cfg, tfr_plcmb_low{:});
tfr_plcmb_high = ft_appendfreq(cfg, tfr_plcmb_high{:});

cfg=[];
cfg.avgoverrpt = 'yes';
cfg.trials = tfr_plcmb_low.trialinfo(:,1)==1;
L_low = ft_selectdata(cfg, tfr_plcmb_low);
L_high = ft_selectdata(cfg, tfr_plcmb_high);
cfg.trials = tfr_plcmb_low.trialinfo(:,1)==2;
R_low = ft_selectdata(cfg, tfr_plcmb_low);
R_high = ft_selectdata(cfg, tfr_plcmb_high);

cfg=[];
cfg.operation = '(x1-x2)./(x1+x2)';
cfg.parameter = 'powspctrm';
ami_low = ft_math(cfg, L_low, R_low);
ami_high = ft_math(cfg, L_high, R_high);


cfg=[];
cfg.avgoverrpt = 'yes';
tfr_plcmb_low = ft_selectdata(cfg, tfr_plcmb_low);
tfr_plcmb_high = ft_selectdata(cfg, tfr_plcmb_high);


if dosave
  save([projectdir, sprintf('results/freq/sub%02d_tfr.mat', subj)], 'tfr_plcmb_low', 'tfr_plcmb_high', 'ami_low', 'ami_high', 'R_low', 'L_low');
end


