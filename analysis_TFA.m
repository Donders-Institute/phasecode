function analysis_TFA(subj,ses)

%   load(sprintf('/home/electromag/matves/Data/phasecode/meg_clean/subject%02d/alldata.mat', subj), 'alldata');
  load(sprintf('/home/electromag/matves/Data/phasecode/meg_clean/subject%02d/session%d.mat', subj, ses), 'cleandata');

% not possible to do time-frequency analysis on the appended data because
% the grad fields are not specified (they cannot be averaged). Better to do
% TFA on the individual sessions and use a grand average.
lag=32.5/1000;


%%%%%%%%%%%%%%%%%%
% Correct trials %
%%%%%%%%%%%%%%%%%%
cleandatacorrect = cleandata;
idx = cleandata.trialinfo(:,7)==1;
cleandatacorrect.trialinfo = cleandatacorrect.trialinfo(idx,:);
cleandatacorrect.trial = cleandatacorrect.trial(:,idx);
cleandatacorrect.time = cleandatacorrect.time(:,idx);



%%%%%%%%%%%%%%%%%%%%
%% Planar gradient %
%%%%%%%%%%%%%%%%%%%%
cfg                 = [];
cfg.method          = 'template';
cfg.template        = 'CTF275_neighb.mat';
neighbours          = ft_prepare_neighbours(cfg, cleandatacorrect);
% neighbours          = ft_prepare_neighbours(cfg, alldatacorrect);
% neighbours          = ft_prepare_neighbours(cfg, alldata);

cfg                 = [];
cfg.method          = 'sincos';
cfg.neighbours      = neighbours;
data_planar = ft_megplanar(cfg, cleandatacorrect);
% data_planar = ft_megplanar(cfg, alldatacorrect);
% data_planar         = ft_megplanar(cfg, alldata);


%% TFR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TFR, Hanning window, fixed window length %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for low frequencies (<30 Hz)
cfg              = [];
cfg.output       = 'pow';
% cfg.channel      = {'MZO' 'MLO' 'MRO'};
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 2:2:30;                         % analysis 2 to 30 Hz in steps of 2 Hz 
cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
cfg.toi          = -0.8+lag:0.05:1.25+lag;                  % time window "slides" from -0.8 to 1.25 sec in steps of 0.05 sec (50 ms)
cfg.keeptrials = 'yes';
TFRhann = ft_freqanalysis(cfg, data_planar);
% TFRhann = ft_freqanalysis(cfg, alldatacorrect);
cfg=[];
tl_plancmbLF          = ft_combineplanar(cfg, TFRhann);



%%%%%%%%%%%%%%%%%%%
% TFR, Multitaper %
%%%%%%%%%%%%%%%%%%%
% for high (30<) frequencies)
cfg = [];
cfg.output     = 'pow';
% cfg.channel    = {'MZO' 'MLO' 'MRO'}; %{'MLO' 'MRO'};
cfg.method     = 'mtmconvol';
cfg.foi        = 30:2:80;
cfg.t_ftimwin  = 5./cfg.foi;
cfg.tapsmofrq  = 0.4 *cfg.foi; %width of frequency smoothing
cfg.toi        = -0.8+lag:0.05:1.25+lag;
cfg.keeptrials = 'yes';
TFRmult = ft_freqanalysis(cfg, data_planar);
% TFRmult = ft_freqanalysis(cfg, alldatacorrect);
cfg=[];
tl_plancmbHF          = ft_combineplanar([], TFRmult);






%% Phase estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TFR, Hanning window, fixed window length %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfg = [];
cfg.method = 'mtmconvol';
cfg.keeptrials = 'yes';
cfg.taper = 'hanning'; 
cfg.toi = -0.8+lag:0.05:1.25+lag; % assuming your trials are of equal length
cfg.foi = 8:12; % or whatever frequency you want estimated (can of course be a vector)
cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5; % use e.g. 5 cycles to estimate the phase
cfg.output = 'fourier';
Fourier = ft_freqanalysis(cfg, cleandatacorrect);

Fourier.pow = abs(Fourier.fourierspctrm);
Fourier.phs = angle(Fourier.fourierspctrm);

save(sprintf('TFA_%d_%d', subj, ses), 'tl_plancmbLF', 'tl_plancmbHF', 'Fourier', '-v7.3');
end

