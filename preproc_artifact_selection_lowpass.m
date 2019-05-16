function preproc_artifact_selection_lowpass(subj, ses)
datainfo; %load subject info

dt=1/400;
beamerLag = (2/60);
lag = dt*round(beamerLag/dt);

%% define trials
cfg = subjects(subj);
cfg.dataset = subjects(subj).session(ses).dataset;
cfg.trialfun = subjects(subj).trialfun;
cfg = ft_definetrial(cfg);
format('shortG') %make the cfg.trl more readible

% read and preprocess data (filtering, baseline correction etc)
cfg.channel = subjects(subj).channels;
cfg.continuous = 'yes';
cfg.bpfilter = 'yes';
cfg.bpfreq = [0.1 40];
cfg.bpfiltord = 2;
data = ft_preprocessing(cfg);

%
%% Load artifacts
% cfg=[];
load(sprintf('/home/electromag/matves/Data/phasecode/artifacts/subject%02d/session%d.mat', subj, ses));
% cfg.artfctdef.eye.artifact = artfctdef.eye.artifact;
% cfg.artfctdef.badtrial = artfctdef.badtrial;

%% Reject artifacts from complete data
cfg = [];
cfg.artfctdef = artfctdef;
cfg.artfctdef.reject = 'complete'; % remove complete trials
cfg.artfctdef.crittoilim = [-0.3+lag 1+lag]; % only remove trials that have an artifact within first second (grating) or just before the grating onset
dataclean = ft_rejectartifact(cfg, data);

%% Downsample
cfg=[];
cfg.detrend = 'no';
cfg.demean = 'no';
cfg.resamplefs = 200;
datacleanresampled = ft_resampledata(cfg,dataclean);

save(sprintf('/home/electromag/matves/Data/phasecode/meg_clean/subject%02d/session%ddownsample200_lowpass.mat', subj, ses), 'datacleanresampled');
%}

%%%%%%%
% ICA %
%%%%%%%

%{
%% definetrial for all MEG sensors
cfg = subjects(subj);
cfg.dataset = subjects(subj).session(ses).dataset;
cfg.trialfun = subjects(subj).trialfun;
cfg = ft_definetrial(cfg);
format('shortG') %make the cfg.trl more readible

% read and preprocess data (filtering, baseline correction etc)
cfg.channel = 'MEG';
cfg.continuous = 'yes';
cfg.bpfilter = 'yes';
cfg.bpfreq = [0.1 40];
cfg.bpfiltord = 2;
dataMEG = ft_preprocessing(cfg);

%% ICA comp
cfg = [];
cfg.resamplefs = 140;
cfg.detrend    = 'no';
datads = ft_resampledata(cfg, dataMEG);
% perform the independent component analysis (i.e., decompose the data)
cfg        = [];
cfg.method = 'fastica';
cfg.channel = 'MEG';
cfg.fastica.numOfIC = 50;
comp = ft_componentanalysis(cfg, datads);
save(sprintf('/home/electromag/matves/Data/phasecode/artifacts/subject%02d/session%dcomp_lowpass.mat', subj, ses), 'comp');
%}

% qsubfeval(@preproc_artifact_rejection_lowpass, subj,ses, 'timreq', 3600, 'memreq', 20*1024^3)


end