ft_hastoolbox('dmlt',1);
[v,ftpath] = ft_version;
addpath(fullfile(ftpath,'external/dmlt/external/glmnet'));
projectdir = '/project/301105.02/phasecode';

hemifield = 1;
class = 1;
tBegin = 0;
subj=2;
Fs=200;

datainfo;
for ses=1:3
  filename = [projectdir, sprintf('data/sub%02d/sub%02d-meg%02d/sub%02d-meg%02d_cleandata.mat', subj, subj, ses, subj, ses)];
  tmp{ses} = load(filename, 'data');
  tmp{ses} = tmp{ses}.data;
end
alldata = ft_appenddata([], tmp{:});

% select channels and trials first, to save processing time later on
idx_val_trial = find((alldata.trialinfo(:,2)==14) & alldata.trialinfo(:,7)==1 |...
  (alldata.trialinfo(:,2)==11) & (alldata.trialinfo(:,7)==1));
cfg        = [];
cfg.trials = idx_val_trial;
cfg.channel= {'MLO', 'MRO','MZO'};
alldata = ft_selectdata(cfg, alldata);

cfg            =[];
cfg.lpfilter   = 'yes';
cfg.lpfreq     = 40;
cfg.lpfilttype = 'firws';
cfg.usefftfilt = 'yes';
alldata        = ft_preprocessing(cfg, alldata);


%Hmmm, some hocus pocus I don't understand 
dt   =1/Fs;
lag  = 35/1000;
time = (-0.3+dt):dt:(1.2-0.02+dt);
time = time + lag;
for trl=1:length(alldata.time)
  alldata.time{trl} = alldata.time{trl}+0.0025; % why is lag 35 ms and here 25 ms added?
end

validdata = alldata;

trl_idxCCW = (validdata.trialinfo(:,2)==14); %& alldata.trialinfo(:,7)==1;
trl_idxCW =  (validdata.trialinfo(:,2)==11);% & (alldata.trialinfo(:,7)==1);
cor_val_trlCCW = find(trl_idxCCW)'; % all correct and validly cued trial indices (left cued)
cor_val_trlCW = find(trl_idxCW)';

cfg         = [];
cfg.latency = [-0.1 1.20];
cfg.trials  = cor_val_trlCCW;
dataCCW     = ft_selectdata(cfg, validdata);
cfg.trials  = cor_val_trlCW;
dataCW      = ft_selectdata(cfg, validdata);
    
ntrialCW  = length(dataCW.trial);
ntrialCCW = length(dataCCW.trial);
ntrials   = min(ntrialCW, ntrialCCW);

cfg        = [];
cfg.trials = 1:ntrials;
dataCW     = ft_selectdata(cfg,dataCW);
dataCCW    = ft_selectdata(cfg,dataCCW);
nchan      = min(size(dataCW.trial{1}));

time = dataCW.time{1};


%% Train
% now for every window of interest (WOI) train on all timepoints of all trials (except TOI) of all trials and test on window of interest from 1 trial.
% for this, put all windows of all trials underneath each other
idx = find(roundn(time,-4)==(roundn(tBegin,-4)));
modelTested = zeros(522,2,ntrials);

for trial = 1:ntrials %-> hard coded pairs the CW and CCW test trials, is this intended, or does it matter
  trial
    cfg       = [];
  cfg.trials = trial;
  testCW    = ft_selectdata(cfg, dataCW);
  testCCW   = ft_selectdata(cfg, dataCCW);
  
  cfg         = [];
  cfg.trials  = setdiff(1:numel(dataCW.trial),trial);
  cfg.latency = [0.2 inf];
  trainCW    = ft_selectdata(cfg, dataCW); 
  trainCCW   = ft_selectdata(cfg, dataCCW);
  
  nTrainWindow    = 1;%5;
  cfg = [];
  cfg.length = nTrainWindow./Fs;
  trainCW   = ft_redefinetrial(cfg, trainCW);
  trainCCW  = ft_redefinetrial(cfg, trainCCW);
%   cfg.overlap = 0.8;
  testCW  = ft_redefinetrial(cfg, testCW);
  testCCW = ft_redefinetrial(cfg, testCCW);
  
  
  nTrain = numel(trainCW.trial);%3000;
  cfg = [];
  cfg.trials = randperm(numel(trainCW.trial),nTrain);
  trainCWsub   = ft_selectdata(cfg, trainCW);
  cfg.trials = randperm(numel(trainCCW.trial),nTrain);
  trainCCWsub  = ft_selectdata(cfg, trainCCW);
  
  TrainData = cat(3, trainCWsub.trial{:}, trainCCWsub.trial{:});
  TrainData = reshape(TrainData,[], nTrain*2)';

  % Z-scoring (mean and std of only training data, not test data)
  mean_data = mean(TrainData);
  std_data  = std(TrainData,[],1);
  TrainData = (TrainData-mean_data)./std_data;
  TrainDesign = [ones(nTrain,1); 2*ones(nTrain,1)];
  %model = dml.enet('family', 'binomial', 'df', 0, 'alpha', 0.01, 'nlambda', 50 );
  %model = model.train(TrainData,TrainDesign); 
  [B,dev,stats]=mnrfit(TrainData,TrainDesign);
  
  %% Test
  TestDataCW = cat(3, testCW.trial{:});
  TestDataCW = reshape(TestDataCW,[],numel(testCW.trial))';
  
  TestDataCCW = cat(3, testCCW.trial{:});
  TestDataCCW = reshape(TestDataCCW,[],numel(testCCW.trial))';
  
%   TestData = cat(1, TestDataCW, TestDataCCW);
%   TestData = (TestData-mean(TestData))./std(TestData,[],1);
  
  TestDataCW  = (TestDataCW-mean_data)./std_data;
  TestDataCCW = (TestDataCCW-mean_data)./std_data;
  TestData = [TestDataCW; TestDataCCW];
  yhat=mnrval(B,TestData);
  
  modelTested(:,:,trial)= yhat;
end

% sprintf('preprocessing: %d s, processing in loop: %d s per trial', t2, t3)
% sprintf('training model: %d s, testing model: %d', t4, t5+t6)
