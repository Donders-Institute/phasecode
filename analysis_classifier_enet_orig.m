function analysis_classifier_enet_orig(subj, tBegin, class, hemifield)

% class=1 means orientation of both gratings (CW_CW vs CCW_CCW)
% class=2 means orientation of attended grating (CW_CW vs CCW_CCW)
% class=3 means attend left vs attend right irrespective of orientations
if nargin<3
    class = 2;
end
if isempty(class)
    class = 2;
end

if nargin<4
    hemifield = 1; % left hemifield
end
if isempty(hemifield)
    hemifield = 1;
end
Fs=200;
% addpath(genpath('/home/common/matlab/fieldtrip/external/dmlt')); % make sure the multivariate toolbox is in your path
% load(sprintf('/home/electromag/matves/Data/phasecode/meg_clean/subject%02d/alldata%d_lowpass_40.mat', subj, Fs), 'alldata');
datainfo;
for k=1:numel(subjects(subj).sessions)
    ses = subjects(subj).sessions(k);
%     filename = [datadir, sprintf('3016045.07_matves_%03d_%03d', subj, ses), '/cleandata.mat'];
    filename = [projectdir, sprintf('data/sub%02d/sub%02d-meg%02d/sub%02d-meg%02d_cleandata.mat', subj, subj, ses, subj, ses)];
    tmp{k} = load(filename, 'data');
    tmp{k} = tmp{k}.data;
end
alldata = ft_appenddata([], tmp{:});

% cfg=[];
% cfg.lpfilter = 'yes';
% cfg.lpfreq = 40;
% cfg.lpfilttype = 'firws';
% alldata = ft_preprocessing(cfg, alldata);

% dt = 1/400; % 1/Fs
dt=1/Fs;
% lag=35/1000;
if Fs==50 || Fs==100
    lag=40/1000;
elseif Fs==200
    lag=35/1000;
elseif Fs==150 || Fs==120
    lag=2/60;
end
% Use a similar procedure to Ramkumar et al (2013): Use a moving 20ms
% window and slide it over every 2.5ms.
% beamerLag = (2/60); % 2 frames delay, 60Hz refresh (-->33.3 ms)
% lag = 35/1000; % In frames of 5ms-->35ms
time = (-0.3+dt):dt:(1.2-0.02+dt); % start time resolved classification 300ms prestim, end 200 ms poststim (260 in Ramkumar et al)
% the last time window starts at 1.1825 seconds, so that the end of the
% last window is at 1.2s
time = time + lag; % account for 2 frame beamerlag


%%
%%%%%%%%%%%%%%%%%%%
% Trial selection %
%%%%%%%%%%%%%%%%%%%
for trl=1:length(alldata.time)
    alldata.time{trl} = alldata.time{trl}+0.0025;
end

if class==1 %
    % select all trials with CCW_CCW or CW_CW orientation
    idx_val_trial = find((alldata.trialinfo(:,2)==14) & alldata.trialinfo(:,7)==1 |...
        (alldata.trialinfo(:,2)==11) & (alldata.trialinfo(:,7)==1));
    cfg=[];
    cfg.trials = idx_val_trial;
    cfg.channel = 'MEG';%{'MLO', 'MRO','MZO'};
    validdata = ft_selectdata(cfg,alldata);
    
    trl_idxCCW = validdata.trialinfo(:,2)==14 & validdata.trialinfo(:,7)==1 & validdata.trialinfo(:,1)==validdata.trialinfo(:,4);
    trl_idxCW =  validdata.trialinfo(:,2)==11 & (validdata.trialinfo(:,7)==1) & validdata.trialinfo(:,1)==validdata.trialinfo(:,4);
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
elseif class==2
    
    % select all valid cued correct trials. Seperate them for left and right
    % cued
    if hemifield==1 % left
        trl_idxL = (alldata.trialinfo(:,1)==2) & alldata.trialinfo(:,7)==1 & alldata.trialinfo(:,4)==2;
        cor_val_trlL = find(trl_idxL)'; % all correct and validly cued trial indices (left cued)
        cfg=[];
        cfg.channel = {'MLO', 'MRO','MZO'};
        cfg.trials = cor_val_trlL;
        validdataL = ft_selectdata(cfg,alldata);
        
        
        % now seperate them on orientation as well.
        trl_idxL_CCW = find(validdataL.trialinfo(:,2)==13 | validdataL.trialinfo(:,2)==14);
        trl_idxL_CW = find(validdataL.trialinfo(:,2)==11 | validdataL.trialinfo(:,2)==12); % cued grating (left) is CW
        
        cfg = [];
        %         cfg.latency = [0+lag 1+lag]; %[-0.3+lag 1.26+lag];
        cfg.latency = [-0.3+dt+lag 1.20+lag];
        cfg.trials = trl_idxL_CW;
        dataCW = ft_selectdata(cfg, validdataL);
        cfg.trials = trl_idxL_CCW;
        dataCCW = ft_selectdata(cfg, validdataL);
        
        ntrialCW = length(dataCW.trial);
        ntrialCCW = length(dataCCW.trial);
        ntrials = min(ntrialCW, ntrialCCW);
        
        cfg=[];
        cfg.trials = 1:ntrials;
        dataCW = ft_selectdata(cfg,dataCW);
        dataCCW = ft_selectdata(cfg,dataCCW);
        nchan = min(size(dataCW.trial{1}));
    elseif hemifield==2 %right
        trl_idxR =     (alldata.trialinfo(:,1)==1) & (alldata.trialinfo(:,7)==1) & (alldata.trialinfo(:,4)==1);
        cor_val_trlR = find(trl_idxR)';
        cfg=[];
        cfg.trials = cor_val_trlR;
        cfg.channel = {'MLO', 'MRO','MZO'};
        validdataR = ft_selectdata(cfg,alldata);
        
        % now seperate them on orientation as well.
        trl_idxR_CCW = find(validdataR.trialinfo(:,2)==12 | validdataR.trialinfo(:,2)==14); % cued grating (right) is CCW
        trl_idxR_CW = find(validdataR.trialinfo(:,2)==11 | validdataR.trialinfo(:,2)==13); % cued grating (left) is CW
        
        cfg = [];
        %         cfg.latency = [0+lag 1+lag]; %[-0.3+lag 1.26+lag];
        cfg.latency = [-0.3+dt+lag 1.20+lag];
        cfg.trials = trl_idxR_CW;
        dataCW = ft_selectdata(cfg, validdataR);
        cfg.trials = trl_idxR_CCW;
        dataCCW = ft_selectdata(cfg, validdataR);
        
        
        ntrialCW = length(dataCW.trial);
        ntrialCCW = length(dataCCW.trial);
        ntrials = min(ntrialCW, ntrialCCW);
        
        cfg=[];
        cfg.trials = 1:ntrials;
        dataCW = ft_selectdata(cfg,dataCW);
        dataCCW = ft_selectdata(cfg,dataCCW);
        nchan = min(size(dataCW.trial{1}));
    end
    
elseif class==3 %
    % select all correct trials with CCW_CCW or CW_CW orientation
    idx_val_trial = alldata.trialinfo(:,7)==1 & alldata.trialinfo(:,1)==alldata.trialinfo(:,4);
    cfg=[];
    cfg.trials = idx_val_trial;
    cfg.channel = {'MLO', 'MRO','MZO'};
    validdata = ft_selectdata(cfg,alldata);
    
    
    trl_idxR = (validdata.trialinfo(:,4)==2); %& cue right!
    trl_idxL =  (validdata.trialinfo(:,4)==1);% cue left!
    cor_val_trlR = find(trl_idxR)'; % all correct and validly cued trial indices (left cued)
    cor_val_trlL = find(trl_idxL)';
    
    cfg=[];
    cfg.latency = [-0.3+dt+lag 1.20+lag];
    cfg.trials = cor_val_trlR;
    dataCCW = ft_selectdata(cfg, validdata); %actually not CCW/CW but Left/Right!!!!
    cfg.trials = cor_val_trlL;
    dataCW = ft_selectdata(cfg, validdata);
    
    ntrialCW = length(dataCW.trial);
    ntrialCCW = length(dataCCW.trial);
    ntrials = min(ntrialCW, ntrialCCW);
    
    cfg=[];
    cfg.trials = 1:ntrials;
    dataCW = ft_selectdata(cfg,dataCW);
    dataCCW = ft_selectdata(cfg,dataCCW);
    nchan = min(size(dataCW.trial{1}));
end


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
for trial = 1:ntrials
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
    mean_data = repmat(mean(TrainData,1), 2*(ntrials-1),1);
    std_data = repmat(std(TrainData, [], 1), 2*(ntrials-1),1);
    TrainData = (TrainData-mean_data)./std_data;
    TrainDesign = [ones(size(tmpTrainDataCW,1),1); 2*ones(size(tmpTrainDataCCW,1),1)];
    model = dml.enet('family', 'binomial', 'df', 0, 'alpha', 0.1 );
    % model = dml.analysis({dml.svm});%dml.analysis({dml.standardizer dml.svm}); % use svm classification, you can also use enet by putting dml.enet
    % model = model.train(TrainData([1:50000 200001:250000],:),TrainDesign([1:50000 200001:250000], :)); % here you create the model with which you are going to classify your test trials
%     keyboard
    model = model.train(TrainData,TrainDesign); % here you create the model with which you are going to classify your test trials
    
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
    TestDataCW = (TestDataCW-mean_data(1,:))./std_data(1,:);
    TestDataCCW = (TestDataCCW-mean_data(1,:))./std_data(1,:);
    
    modelTested(trial,:)= model.test(TestDataCW);
    modelTested(trial+ntrials,:)= model.test(TestDataCCW);
end
%% Save
% if class==1
%     save(sprintf('/home/electromag/matves/Results/ENET/subject%02d/enet_trial_window_%d_%d_MZO.mat', subj, Fs, idx),'modelTested', 'ntrials');
% elseif class==2
%     if hemifield==1
%         save(sprintf('/home/electromag/matves/Results/ENET/subject%02d/enet_trial_windowL_%d_%d_MZO.mat', subj,Fs, idx),'modelTested', 'ntrials');
%     elseif hemifield==2
%         save(sprintf('/homei/electromag/matves/Results/ENET/subject%02d/enet_trial_windowR_%d_%d_MZO.mat', subj,Fs, idx),'modelTested', 'ntrials');
%     end
% elseif class==3
%     save(sprintf('/home/electromag/matves/Results/ENET/subject%02d/enet_trialCue_window_%d_%d_MZO.mat', subj, Fs, idx),'modelTested','ntrials');
% end
filename = sprintf('/project/3011085.02/phasecode/results/enet/sub%02d/sub%02d_enet_orig_%d.mat', subj, subj, idx);
save(filename, 'modelTested', 'ntrials')

end