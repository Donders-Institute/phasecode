function preproc_artifact_selection(subj, ses)
datainfo; %load subject info
subj=10;ses=22;
dt=1/400;
beamerLag = (2/60);
lag = 32.5/1000;
%
%% define trials
cfg = subjects(subj);
cfg.dataset = subjects(subj).session(ses).dataset;
cfg.trialfun = subjects(subj).trialfun;
cfg = ft_definetrial(cfg);
format('shortG') %make the cfg.trl more readible

% read and preprocess data (filtering, baseline correction etc)
cfg.channel = subjects(subj).channels;
cfg.continuous = 'yes';
cfg.demean = 'yes'; %subtract mean from each channel
cfg.bsfilter = 'yes';
cfg.bsfreq = [49 51; 99 101; 149 151];
data = ft_preprocessing(cfg);

%% Show summary
%
cfg          = [];
cfg.method   = 'summary';
cfg.channel = {'MEG', 'UADC005', 'UADC006'};
cfg.layout = 'CTF275.lay'; %so timecourses and topographies of individual trials can be looked at
% cfg.alim     = 1e-12; 
dataRejVis        = ft_rejectvisual(cfg,data); 
% you can check the rejected trial number by typing
for i=1:length(dataRejVis.cfg.artfctdef.summary.artifact)
trlind(i)=find(data.sampleinfo(:,1)==dataRejVis.cfg.artfctdef.summary.artifact(i));
end;
disp(trlind);
data.badtrials=trlind;
% save data data %could be used to save the bad trialnumbers with the raw
% data.
%

%% Find jump artifacts - not necessary anymore
%{
cfg = [];
 
% channel selection, cutoff and padding
cfg.artfctdef.zvalue.channel    = 'MEG';
cfg.artfctdef.zvalue.cutoff     = 20; % artifact threshold
cfg.artfctdef.zvalue.trlpadding = 0; 
cfg.artfctdef.zvalue.artpadding = 0;
cfg.artfctdef.zvalue.fltpadding = 0;
 
% algorithmic parameters
cfg.artfctdef.zvalue.cumulative    = 'yes';
cfg.artfctdef.zvalue.medianfilter  = 'yes';
cfg.artfctdef.zvalue.medianfiltord = 9;
cfg.artfctdef.zvalue.absdiff       = 'yes';
 
% make the process interactive
cfg.artfctdef.zvalue.interactive = 'yes';
 
[~, artifact_jump] = ft_artifact_zvalue(cfg, data);
%}

%% Find muscle artifacts - not necessary anymore
%{
cfg            = [];

% channel selection, cutoff and padding
cfg.artfctdef.zvalue.channel = 'MEG';
cfg.artfctdef.zvalue.cutoff      = 4; % artifact threshold
cfg.artfctdef.zvalue.trlpadding  = 0; % pad at both sides of trial till it is ... seconds.
cfg.artfctdef.zvalue.fltpadding  = 0;
cfg.artfctdef.zvalue.artpadding  = 0.2; % take 0.2s padding of an artifact (of the part that is above threshold)

% algorithmic parameters
cfg.artfctdef.zvalue.bpfilter    = 'yes';
cfg.artfctdef.zvalue.bpfreq      = [110 140];
cfg.artfctdef.zvalue.bpfiltord   = 9;
cfg.artfctdef.zvalue.bpfilttype  = 'but';
cfg.artfctdef.zvalue.hilbert     = 'yes';
cfg.artfctdef.zvalue.boxcar      = 0.2;

% make the process interactive
cfg.artfctdef.zvalue.interactive = 'yes';

[~, artifact_muscle] = ft_artifact_zvalue(cfg, data);
%}

%% Find EOG artifacts
%{
subj=7;ses=1;
% first define trials only in critical period
cfg            = [];
cfg.dataset = subjects(subj).session(ses).dataset;
cfg.trialfun = 'mytrialfun_critical_period';
cfg = ft_definetrial(cfg);
cfg.datafile = subjects(subj).session(ses).dataset;
cfg.headerfile = subjects(subj).session(ses).dataset;
cfg.continuous = 'yes';

% channel selection, cutoff and padding
cfg.artfctdef.zvalue.channel     = {'UADC005','UADC006'}; % Eye tracker data
% cfg.artfctdef.zvalue.channel     = {'EEG...','EEG...'}; % EOG data.
% Should be one of EEG057 - EEG064
cfg.artfctdef.zvalue.cutoff      = 2;
cfg.artfctdef.zvalue.trlpadding  = 0;
cfg.artfctdef.zvalue.artpadding  = 0.1;
cfg.artfctdef.zvalue.fltpadding  = 0;

% algorithmic parameters
cfg.artfctdef.zvalue.bpfilter   = 'yes';
cfg.artfctdef.zvalue.bpfilttype = 'but';
cfg.artfctdef.zvalue.bpfreq     = [1 15];
cfg.artfctdef.zvalue.bpfiltord  = 4;
cfg.artfctdef.zvalue.hilbert    = 'yes';

% feedback
cfg.artfctdef.zvalue.interactive = 'yes';

[~, artifact_EOG] = ft_artifact_zvalue(cfg);
artifact_visual = artifact_EOG;
% Save artifacts
load(sprintf('/home/electromag/matves/Data/phasecode/artifacts/subject%02d/session%d.mat', subj, ses));
cfg = [];
cfg.channel     = 'MEG';
% cfg.artfctdef.eye.artifact = artifact_EOG;

% cfg.artfctdef.jump.artifact = artifact_jump;
% cfg.artfctdef.muscle.artifact = artifact_muscle;
% cfg.artfctdef.badtrial = trlind;
artfctdef.visual.artifact = artifact_visual; %these are the artifacts in the prestim period
save(sprintf('/home/electromag/matves/Data/phasecode/artifacts/subject%02d/session%d.mat', subj, ses), 'artfctdef');
%}

%% when loading previously defined artifacts:
% cfg=[];
% % cfg.channel=[];
load(sprintf('/home/electromag/matves/Data/phasecode/artifacts/subject%02d/session%d.mat', subj, ses));
% cfg.artfctdef.eye.artifact = artfctdef.eye.artifact;
% cfg.artfctdef.badtrial = artfctdef.badtrial;

% Reject artifacts from complete data
cfg = [];
cfg.artfctdef = artfctdef;
cfg.artfctdef.reject = 'complete'; % remove complete trials
cfg.artfctdef.crittoilim = [-0.3+lag 1+lag]; % only remove trials that have an artifact within first second (grating) or just before the grating onset
dataclean = ft_rejectartifact(cfg, tmp);

%% Downsample
cfg=[];
cfg.detrend = 'no';
cfg.demean = 'no';
cfg.resamplefs = 400;
datacleanresampled = ft_resampledata(cfg,dataclean);

% save(sprintf('/home/electromag/matves/Data/phasecode/meg_clean/subject%02d/session%ddownsample.mat', subj, ses), 'datacleanresampled');
%{

%% define trials for all MEG sessions
cfg = subjects(subj);
cfg.dataset = subjects(subj).session(ses).dataset;
cfg.trialfun = subjects(subj).trialfun;
cfg = ft_definetrial(cfg);
format('shortG') %make the cfg.trl more readible

% read and preprocess data (filtering, baseline correction etc)
cfg.channel = 'MEG';
cfg.continuous = 'yes';
cfg.bsfilter = 'yes';
cfg.bsfreq = [49 51; 99 101; 149 151];
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
save(sprintf('/home/electromag/matves/Data/phasecode/artifacts/subject%02d/session%dcomp.mat', subj, ses), 'comp');
%}
% qsubfeval(@preproc_artifact_rejection, subj, ses, 'timreq', 3600, 'memreq', 20*1024^3)
end
