    addpath('/project/3011085.02/scripts/fieldtrip/external/dmlt/')
addpath('/project/3011085.02/scripts/fieldtrip/external/dmlt/external/')
    hemifield = 1;
    class = 1;
    tBegin = 0;
    subj=4;
Fs=200;
datainfo;
t1 = tic;
for k=1:numel(subjects(subj).sessions)
    ses = subjects(subj).sessions(k);
    filename = [projectdir, sprintf('data/sub%02d/sub%02d-meg%02d/sub%02d-meg%02d_cleandata.mat', subj, subj, ses, subj, ses)];
    tmp{k} = load(filename, 'data');
    tmp{k} = tmp{k}.data;
end
alldata = ft_appenddata([], tmp{:});

cfg=[];
cfg.lpfilter = 'yes';
cfg.lpfreq = 40;
cfg.lpfilttype = 'firws';
alldata = ft_preprocessing(cfg, alldata);

dt=1/Fs;
lag=35/1000;


time = (-0.3+dt):dt:(1.2-0.02+dt);
time = time + lag;


%%
%%%%%%%%%%%%%%%%%%%
% Trial selection %
%%%%%%%%%%%%%%%%%%%
for trl=1:length(alldata.time)
    alldata.time{trl} = alldata.time{trl}+0.0025;
end

    % select all trials with CCW_CCW or CW_CW orientation
    idx_val_trial = find((alldata.trialinfo(:,2)==14) & alldata.trialinfo(:,7)==1 |...
        (alldata.trialinfo(:,2)==11) & (alldata.trialinfo(:,7)==1));
    cfg=[];
    cfg.trials = idx_val_trial;
    cfg.channel = {'MLO', 'MRO','MZO'};
    validdata = ft_selectdata(cfg,alldata);
    
    trl_idxCCW = (validdata.trialinfo(:,2)==14); %& alldata.trialinfo(:,7)==1;
    trl_idxCW =  (validdata.trialinfo(:,2)==11);% & (alldata.trialinfo(:,7)==1);
    cor_val_trlCCW = find(trl_idxCCW)'; % all correct and validly cued trial indices (left cued)
    cor_val_trlCW = find(trl_idxCW)';
    
    cfg=[];
    cfg.latency = [-0.3+dt+lag 1.20+lag];
    cfg.trials = cor_val_trlCCW;
    dataCCW = ft_selectdata(cfg, validdata);
    cfg.trials = cor_val_trlCW;
    dataCW = ft_selectdata(cfg, validdata);
    
    ntrialCW = length(dataCW.trial);
    ntrialCCW = length(dataCCW.trial);
    ntrials = min(ntrialCW, ntrialCCW);
    
    cfg=[];
    cfg.trials = 1:ntrials;
    dataCW = ft_selectdata(cfg,dataCW);
    dataCCW = ft_selectdata(cfg,dataCCW);
    nchan = min(size(dataCW.trial{1}));

%%
%%%%%%%%%%%%
% Classify %
%%%%%%%%%%%%

tBegin=tBegin+lag;


%% Train
% now for every window of interest (WOI) train on all timepoints of all trials (except TOI) of all trials and test on window of interest from 1 trial.
% for this, put all windows of all trials underneath each other
idx = find(roundn(time,-4)==(roundn(tBegin,-4)));
modelTested = zeros(2*ntrials,2);
t2 = toc(t1);
for trial = 1%:ntrials
    CW_noTOI = dataCW;
    CW_TOI = dataCW.trial{trial};
    CW_noTOI.trial(trial)=[];
    CCW_noTOI = dataCCW;
    CCW_TOI = dataCCW.trial{trial};
    CCW_noTOI.trial(trial)=[];
    
    
    % put all windows of all trials underneath eachother
    nTrainWindow = 4;
    tmpTrainDataCW = zeros((ntrials-1), nTrainWindow*nchan);
    tmpTrainDataCCW = zeros((ntrials-1), nTrainWindow*nchan);
    ii=1;
    for trl=1:ntrials-1
        jj=0;
        kk=1;
        for win=1:nTrainWindow
            tmpTrainDataCW(ii,kk:kk+nchan-1) = reshape((CW_noTOI.trial{trl}(:,idx+jj)),1,size(CW_noTOI.trial{trl}(:,idx+jj),1)*(size(CW_noTOI.trial{trl}(:,idx+jj),2)));
            tmpTrainDataCCW(ii,kk:kk+nchan-1) = reshape((CCW_noTOI.trial{trl}(:,idx+jj)),1,size(CCW_noTOI.trial{trl}(:,idx+jj),1)*(size(CCW_noTOI.trial{trl}(:,idx+jj),2)));
            kk=kk+nchan;
            jj=jj+1;
        end
        ii=ii+1;
    end
    
    
    
    TrainData = [tmpTrainDataCW; tmpTrainDataCCW];
    
    % Z-scoring (mean and std of only training data, not test data)
    mean_data = mean2(TrainData);
    std_data = std2(TrainData);
    TrainData = (TrainData-mean_data)./std_data;
    TrainDesign = [ones(size(tmpTrainDataCW,1),1); 2*ones(size(tmpTrainDataCCW,1),1)];
    model = dml.enet('family', 'binomial', 'df', 0, 'alpha', 0.1 );

t3 = toc(t1)-t2;
    model = model.train(TrainData,TrainDesign); % here you create the model with which you are going to classify your test trials
t4 = toc(t1)-t2-t3;
   %% Test
    jj=0;
    kk=1;
    for win=1:nTrainWindow
        TestDataCW(1,kk:kk+nchan-1) = reshape((CW_TOI(:,idx+jj)),1,size(CW_TOI(:,idx+jj),1)*(size(CW_TOI(:,idx+jj),2)));
        TestDataCCW(1,kk:kk+nchan-1) = reshape((CCW_TOI(:,idx+jj)),1,size(CCW_TOI(:,idx+jj),1)*(size(CCW_TOI(:,idx+jj),2)));
        kk=kk+nchan;
        jj=jj+1;
    end
    
    
    % TestData = [TestDataCW; TestDataCCW];
    TestDataCW = (TestDataCW-mean_data)./std_data;
    TestDataCCW = (TestDataCCW-mean_data)./std_data;
t5 = toc(t1)-t2-t3-t4;
    modelTested(trial,:)= model.test(TestDataCW);
    modelTested(trial+ntrials,:)= model.test(TestDataCCW);
t6 = toc(t1)-t2-t3-t4-t5;
end

sprintf('preprocessing: %d s, processing in loop: %d s per trial', t2, t3)
sprintf('training model: %d s, testing model: %d', t4, t5+t6)
