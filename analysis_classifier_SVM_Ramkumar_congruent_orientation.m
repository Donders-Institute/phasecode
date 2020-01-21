function analysis_classifier_SVM_Ramkumar_simple(subj)
% classify on both orientations clockwise vs counterclockwise, irrespective
% of cue (ignore the trials with 1 of each orientation)


addpath(genpath('/home/common/matlab/fieldtrip/external/dmlt')); % make sure the multivariate toolbox is in your path
%     addpath('/home/electromag/matves/Downloads/fieldtrip-20150923/external/dmlt');


% load(sprintf('/home/electromag/matves/Data/phasecode/meg_clean/subject%02d/alldata.mat', subj), 'alldata');


load(sprintf('/home/electromag/matves/Data/phasecode/meg_clean/subject%02d/alldata200_lowpass_40.mat', subj), 'alldata');



dt = 1/200; % 1/Fs
% Use a similar procedure to Ramkumar et al (2013): Use a moving 20ms
% window and slide it over every 2.5ms.
% beamerLag = (2/60); % 2 frames delay, 60Hz refresh (-->33.3 ms)
lag = 35/1000; % In frames of 2.5ms-->32.5ms
time = (-0.3+dt):dt:(1.2-0.02+dt); % start time resolved classification 300ms prestim, end 200 ms poststim (260 in Ramkumar et al)
% the last time window starts at 1.1825 seconds, so that the end of the
% last window is at 1.2s
time = time + lag; % account for 2 frame beamerlag
size_of_window = 20/1000; %10 ms
timepoint_per_window = size_of_window/dt;
%%%%%%%%%%%%%%%%%%%%%%
%%% Z-scoring data %%%
%%%%%%%%%%%%%%%%%%%%%%
% for every time-point, subtract the mean over all trials and devide by the
% standard deviation [sqrt(mean(deviation^2))]. To this for every time
% point and channel seperately.

%%
%%%%%%%%%%%%%%%%%%%
% Trial selection %
%%%%%%%%%%%%%%%%%%%
% select all correct trials with CCW_CCW or CW_CW orientation
for trl=1:length(alldata.time)
    alldata.time{trl} = alldata.time{trl}+dt/2; % realign time axis
end


trl_idxCCW = (alldata.trialinfo(:,2)==14) & (alldata.trialinfo(:,7)==1) & (alldata.trialinfo(:,1)==alldata.trialinfo(:,4));
trl_idxCW =  (alldata.trialinfo(:,2)==11) & (alldata.trialinfo(:,7)==1) & (alldata.trialinfo(:,1)==alldata.trialinfo(:,4));
cor_val_trlCCW = find(trl_idxCCW)'; % all correct and validly cued trial indices (left cued)
cor_val_trlCW = find(trl_idxCW)';

cfg=[];
cfg.latency = [-0.3+dt+lag 1.2+lag];
cfg.trials = cor_val_trlCCW;
dataCCW = ft_selectdata(cfg, alldata);
cfg.trials = cor_val_trlCW;
dataCW = ft_selectdata(cfg, alldata);


ntrialCW = length(dataCW.trial);
ntrialCCW = length(dataCCW.trial);
ntrials = min(ntrialCW, ntrialCCW);

cfg=[];
cfg.trials = 1:ntrials;
cfg.channel = {'MRO', 'MZO','MLO'};
dataCW = ft_selectdata(cfg,dataCW);
dataCCW = ft_selectdata(cfg,dataCCW);
nchan = min(size(dataCW.trial{1}));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate means and std %
%%%%%%%%%%%%%%%%%%%%%%%%%%%


meandata = zeros(length(time),1); % # of windows
stddata = meandata;
for window = 1:length(time)
    jj=1;
    comb_trials = zeros(nchan,2*ntrials*timepoint_per_window);
    for trl = 1:ntrials
        comb_trials(:, jj:jj+timepoint_per_window-1) = dataCW.trial{trl}(:,window:window+timepoint_per_window-1); % copy the same window of all trials and put them together in one matrix
        comb_trials(:, jj+timepoint_per_window:jj+2*timepoint_per_window-1) = dataCCW.trial{trl}(:,window:window+timepoint_per_window-1); % copy the same window of all trials and put them together in one matrix
        jj=jj+8;
    end
    meandata(window) = mean(mean(comb_trials)); % calculate the mean of all matrix elements (mean of 1 window over all trials and channels)
    stddata(window) = std2(comb_trials); % calculate the std of all matrix elements
end

%% double CCW vs double CW


stat = struct('accuracy', {}, 'binomial', {}, 'contingency', {}); %TimeresolvedMoving...
acc = zeros(length(time),1);
bin = zeros(length(time),1);
jj=1;
for tBegin = time; %20ms window moving 1/fs (1/200)
    
    %%%%%%%%%%%%%%%%%%%%%%
    % Select time window %
    %%%%%%%%%%%%%%%%%%%%%%
    cfg=[];
    cfg.latency = [tBegin tBegin+(timepoint_per_window-1)*dt];
    dCCW = ft_selectdata(cfg,dataCCW);
    dCW = ft_selectdata(cfg,dataCW);
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % Z-score time window %
    %%%%%%%%%%%%%%%%%%%%%%%
    for trl = 1:ntrials
        dCCW.trial{trl} = ((dCCW.trial{trl} - meandata(jj))./ stddata(jj));
        dCW.trial{trl} = ((dCW.trial{trl} - meandata(jj))./ stddata(jj));
    end

    %%%%%%%%%%%%%%%%%%%%%%
    % Time-lock analysis %
    %%%%%%%%%%%%%%%%%%%%%%
    %compute time-locked average and covariance matrix
    
    cfg             = [];
%     cfg.parameter   = 'trial';
    cfg.keeptrials  = 'yes'; % classifiers operate on individual trials
%     cfg.channel     = {'MRO', 'MLO'}; %{'MRO'}; % occipital channels contralateral to cued hemifield only
%     cfg.removemean ='no';
% cfg.vartrllength=0;
    tCCW   = ft_timelockanalysis(cfg,dCCW);
    tCW   = ft_timelockanalysis(cfg,dCW);
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % Time-lock statistics %
    %%%%%%%%%%%%%%%%%%%%%%%%
    % use the default classification procedure;  a standardization of the data
    % (subtraction of the mean and division by the standard deviation),
    % followed by applying a linear support vector machine
    
%     for tBegin = time; %20ms window moving 1/fs (1/400)
    cfg         = [];
    cfg.layout  = 'CTF275.lay';
    cfg.method  = 'crossvalidate';% use n-fold cross-validation
    cfg.nfolds = 5; % number of cross-validation folds (default =5)
    cfg.design  = [ones(size(tCW.trial,1),1); 2*ones(size(tCCW.trial,1),1)]';
    %specify design matrix: vector with 1's for leftCCW, 2 for leftCW
    cfg.latency     = [tBegin tBegin+(timepoint_per_window-1)*dt]; % 20ms window
    cfg.statistic = {'accuracy' 'binomial' 'contingency'};
    % cfg.mva = {dml.standardizer dml.enet('family', 'binomial','alpha', 0.2)};
    TrMstat = ft_timelockstatistics(cfg,tCW, tCCW);
    
    model{jj} = TrMstat.model;
    stat{jj} = TrMstat.statistic;
    acc(jj,1) = TrMstat.statistic.accuracy;
%     bin(jj,1) = TrMstat.statistic.binomial;%p-value
    jj=jj+1;
end
% var = acc.*(1-bin);


%% Average and save
% AvgAcc = mean(acc(122:522,:));



save(sprintf('/home/electromag/matves/Results/SVM/subject%02d/SVMramkumar_simple2.mat', subj),'acc', 'stat', 'time');

end