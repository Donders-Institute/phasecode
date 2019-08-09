function analysis_classifier_enet(subj, tBegin, class, hemifield, trainfullwindow, leaveoutwoi)

% class=1 means orientation of both gratings (CW_CW vs CCW_CCW)
% class=2 means orientation of attended grating (CW vs CCW)
% class=3 means attend left vs attend right irrespective of orientations
if ~exist('class', 'var'), class=4; end
if ~exist('hemifield', 'var'), hemifield=1; end
if ~exist('tBegin', 'var'), tBegin = []; end
if ~exist('trainfullwindow', 'var'), trainfullwindow = false; end
if ~exist('leaveoutwoi', 'var'), leaveoutwoi = false; end

datainfo;
for ses=1:3
%     filename = [datadir, sprintf('3016045.07_matves_%03d_%03d', subj, ses), '/cleandata.mat'];
    filename = [projectdir, sprintf('data/sub%02d/sub%02d-meg%02d/sub%02d-meg%02d_cleandata.mat', subj, subj, ses, subj, ses)];
    tmp{ses} = load(filename, 'data');
    tmp{ses} = tmp{ses}.data;
end
data = ft_appenddata([], tmp{:});

cfg=[];
cfg.lpfilter = 'yes';
cfg.lpfreq = 40;
cfg.lpfilttype = 'firws';
data = ft_preprocessing(cfg, data);

clear tmp
dt=1/200;

% Use a similar procedure to Ramkumar et al (2013): Use a moving 20ms
% window and slide it over every 2.5ms.
% beamerLag = (2/60); % 2 frames delay, 60Hz refresh (-->33.3 ms)
% lag = 35/1000; % In frames of 5ms-->35ms
time = (0.5+dt):dt:1.2; % start time resolved classification 300ms prestim, end 200 ms poststim (260 in Ramkumar et al)
% the last time window starts at 1.1825 seconds, so that the end of the
% last window is at 1.2s


%%
%%%%%%%%%%%%%%%%%%%
% Trial selection %
%%%%%%%%%%%%%%%%%%%

if class==1 %
    %{
    % select all trials with CCW_CCW or CW_CW orientation
    idx_val_trial = find((data.trialinfo(:,2)==14) & data.trialinfo(:,7)==1 |...
        (data.trialinfo(:,2)==11) & (data.trialinfo(:,7)==1));
    cfg=[];
    cfg.trials = idx_val_trial;
    cfg.channel = {'MLO', 'MRO','MZO'};
    validdata = ft_selectdata(cfg,data);
    
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
elseif class==2
    
    % select all valid cued correct trials. Seperate them for left and right
    % cued
    if hemifield==1 % left
        trl_idxL = (data.trialinfo(:,1)==2) & data.trialinfo(:,7)==1 & data.trialinfo(:,4)==2;
        cor_val_trlL = find(trl_idxL)'; % all correct and validly cued trial indices (left cued)
        cfg=[];
        cfg.channel = {'MLO', 'MRO','MZO'};
        cfg.trials = cor_val_trlL;
        validdataL = ft_selectdata(cfg,data);
        
        
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
        trl_idxR =     (data.trialinfo(:,1)==1) & (data.trialinfo(:,7)==1) & (data.trialinfo(:,4)==1);
        cor_val_trlR = find(trl_idxR)';
        cfg=[];
        cfg.trials = cor_val_trlR;
        cfg.channel = {'MLO', 'MRO','MZO'};
        validdataR = ft_selectdata(cfg,data);
        
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
    idx_val_trial = data.trialinfo(:,7)==1 & data.trialinfo(:,1)==data.trialinfo(:,4);
    cfg=[];
    cfg.trials = idx_val_trial;
    cfg.channel = {'MLO', 'MRO','MZO'};
    validdata = ft_selectdata(cfg,data);
    
    
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
    %}
elseif class==4
        trl_idx = (data.trialinfo(:,1)==hemifield) ;
        trl_idx = find(trl_idx)'; % all correct and validly cued trial indices (left cued)
        cfg=[];
        cfg.channel = {'MLO', 'MRO','MZO'};
        cfg.trials = trl_idx;
        validdata = ft_selectdata(cfg,data);
        
        if hemifield==1 % left
            % now seperate them on orientation as well.
            trl_idx_CCW = find(validdata.trialinfo(:,2)==13 | validdata.trialinfo(:,2)==14);
            trl_idx_CW = find(validdata.trialinfo(:,2)==11 | validdata.trialinfo(:,2)==12); % cued grating (left) is CW
        elseif hemifield==2
            trl_idx_CCW = find(validdata.trialinfo(:,2)==12 | validdata.trialinfo(:,2)==14); % cued grating (right) is CCW
            trl_idx_CW = find(validdata.trialinfo(:,2)==11 | validdata.trialinfo(:,2)==13); % cued grating (left) is CW
        end
        
        cfg = [];
        cfg.latency = [-0.3 1.2];
        cfg.trials = trl_idx_CW;
        dataCW = ft_selectdata(cfg, validdata);
        cfg.trials = trl_idx_CCW;
        dataCCW = ft_selectdata(cfg, validdata);
        
        ntrialCW = length(dataCW.trial);
        ntrialCCW = length(dataCCW.trial);
        ntrials = min(ntrialCW, ntrialCCW);
        
        cfg=[];
        tmp = randperm(ntrialCW);
        cfg.trials = tmp(1:ntrials);
        dataCW = ft_selectdata(cfg,dataCW);
        tmp = randperm(ntrialCCW);
        cfg.trials = tmp(1:ntrials);
        dataCCW = ft_selectdata(cfg,dataCCW);
        
        nchan = min(size(dataCW.trial{1}));
end


%% Classify 
windowSize = 1/200%0.02;
ntpoints = windowSize/dt;
cfg=[];
if isempty(tBegin)
    toilim = [dt 1.2];
else
    toilim = [tBegin+dt tBegin+windowSize];
end
cfg.toilim = toilim;
dataCW_short = ft_redefinetrial(cfg, dataCW);
dataCCW_short = ft_redefinetrial(cfg, dataCCW);

if trainfullwindow
    cfg=[];
    cfg.toilim = [0.5+dt 1.2];
    dataCW_fullwindow = ft_redefinetrial(cfg, dataCW);
    dataCCW_fullwindow = ft_redefinetrial(cfg, dataCCW);
end

% idx = find(roundn(time,-4)==(roundn(tBegin,-4)));
% modelTested = zeros(2*ntrials,2);
idx_woi = find(time==tBegin+dt);% index for the window of interest, relative to time

probability = zeros(2*ntrials,1);
trialinfo = [dataCW_short.trialinfo; dataCCW_short.trialinfo];
% Treat different windows as trials. Leave out the entire window of 
% interest so that all trials of the window of interest can be tested based
% on the same model.
if trainfullwindow && leaveoutwoi
    ntrials_lowoi = ntrials*numel(time)/ntpoints;
    tmptrainCW = permute(cat(3,dataCW_fullwindow.trial{:}), [3,1,2]);
    tmptrainCW(:,:, idx_woi:idx_woi+ntpoints-1) = []; % use this if testing on all trials at the same time (training on all windows except woi)
    tmptrainCCW = permute(cat(3,dataCCW_fullwindow.trial{:}), [3,1,2]);
    tmptrainCCW(:,:, idx_woi:idx_woi+ntpoints-1) = [];
    [s1,s2,s3] = size(tmptrainCW);
    trainCW = zeros(ntrials_lowoi,s2,ntpoints);
    trainCCW = zeros(ntrials_lowoi,s2,ntpoints);
    idx1=1;
    for l=1:s3/ntpoints
        trainCW(idx1:idx1+s1-1,:,:) = tmptrainCW(:, :, 4*l-3:4*l);
        trainCCW(idx1:idx1+s1-1,:,:) = tmptrainCCW(:, :, 4*l-3:4*l);
        idx1=idx1+s1;
    end
    
    [s4,s5,s6] = size(trainCW);
    trainCW = reshape(trainCW, [s4, s5*s6]);
    trainCCW = reshape(trainCCW, [s4, s5*s6]);
    trainData = zscore([trainCW; trainCCW], [], 2); % every trial has mean 0 and std 1
    
    design = [ones(ntrials_lowoi,1); 2*ones(ntrials_lowoi,1)];
    
    % Train model
    model = dml.enet('family', 'binomial', 'df', 0, 'alpha', 0.1 );
    model = model.train(trainData,design); % here you create the model with which you are going to classify your test trials
    
    for k=1:ntrials
        k
        % prepare test data
        testCW = reshape(dataCW_short.trial{k}, 1, []);
        testCCW = reshape(dataCCW_short.trial{k}, 1, []);
        testCW = zscore(testCW,[],2);
        testCCW = zscore(testCCW,[],2);
        
        % Test model
        tmp = model.test(testCW);
        probability(k,:) = tmp(1,1);         % save result
        tmp = model.test(testCCW);
        probability(k+ntrials,:) = tmp(1,2);
    end
    
else
    for k = 1:ntrials
        k
        % prepare data
        % Treat different windows as trials. For trial k, only leave out
        % the window of interest. Loops over trials. Very inefficient.
        if trainfullwindow
            ntrials_fullwindow = ntrials*numel(time)/ntpoints-1;
           
            tmptrainCW = permute(cat(3,dataCW_fullwindow.trial{:}), [3,1,2]);
            tmptrainCCW = permute(cat(3,dataCCW_fullwindow.trial{:}), [3,1,2]);
            [s1,s2,s3] = size(tmptrainCW);
            trainCW = zeros(ntrials_fullwindow,s2,ntpoints);
            trainCCW = zeros(ntrials_fullwindow,s2,ntpoints);
            idx1=1;
            for l=1:s3/ntpoints
                if l==idx_woi
                    % skip the window of the trial that will be tested (don't
                    % put it in training data).
                    tmp1 = tmptrainCW(:, :, 4*l-3:4*l);
                    tmp1(k,:,:) = [];
                    trainCW(idx1:idx1+s1-2,:,:) = tmp1;
                    tmp2 = tmptrainCCW(:, :, 4*l-3:4*l);
                    tmp2(k,:,:) = [];
                    trainCCW(idx1:idx1+s1-2,:,:) = tmp2;
                    idx1 = idx1+s1-1;
                else
                    trainCW(idx1:idx1+s1-1,:,:) = tmptrainCW(:, :, 4*l-3:4*l);
                    trainCCW(idx1:idx1+s1-1,:,:) = tmptrainCCW(:, :, 4*l-3:4*l);
                    idx1=idx1+s1;
                end
            end
            [s4,s5,s6] = size(trainCW);
            trainCW = reshape(trainCW, [s4, s5*s6]);
            trainCCW = reshape(trainCCW, [s4, s5*s6]);
            trainData = zscore([trainCW; trainCCW], [], 2); % every trial has mean 0 and std 1
            
            design = [ones(ntrials_fullwindow,1); 2*ones(ntrials_fullwindow,1)];
            
        else
            % train the model only on the window of interest of training
            % trials
            trainCW = reshape(permute(cat(3,dataCW_short.trial{:}), [3,1,2]), ntrials, []);
            trainCW(k,:) = [];
            
            
            trainCCW = reshape(permute(cat(3,dataCCW_short.trial{:}), [3,1,2]), ntrials, []);
            trainCCW(k,:) = [];
            
            trainData = zscore([trainCW; trainCCW], [], 2); % every trial has mean 0 and std 1
            design = [ones(ntrials-1,1); 2*ones(ntrials-1,1)];
        end

        % Train model
        model = dml.enet('family', 'binomial', 'df', 0, 'alpha', 0.1 );
        model = model.train(trainData,design); % here you create the model with which you are going to classify your test trials
        
        % Prepare test data
        testCW = reshape(dataCW_short.trial{k}, 1, []);
        testCCW = reshape(dataCCW_short.trial{k}, 1, []);
        testCW = zscore(testCW,[],2);
        testCCW = zscore(testCCW,[],2);
        
        % Test model
        tmp = model.test(testCW);
        tmp(1,1)
        probability(k,:) = tmp(1,1);         % save result
        tmp = model.test(testCCW);
        probability(k+ntrials,:) = tmp(1,2);
        tmp(1,1)
    end
end
 %{   
    CW_noTOI = dataCW;
    CW_TOI = dataCW.trial{k};
    CW_noTOI.trial(k)=[];
    CCW_noTOI = dataCCW;
    CCW_TOI = dataCCW.trial{k};
    CCW_noTOI.trial(k)=[];
    
    
    % put all windows of all trials underneath each other
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
    % model = dml.analysis({dml.svm});%dml.analysis({dml.standardizer dml.svm}); % use svm classification, you can also use enet by putting dml.enet
    % model = model.train(TrainData([1:50000 200001:250000],:),TrainDesign([1:50000 200001:250000], :)); % here you create the model with which you are going to classify your test trials
    keyboard
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
    TestDataCW = (TestDataCW-mean_data)./std_data;
    TestDataCCW = (TestDataCCW-mean_data)./std_data;
    
    modelTested(k,:)= model.test(TestDataCW);
    modelTested(k+ntrials,:)= model.test(TestDataCCW);
end
   %} 
%% Save
filename = sprintf('/project/3011085.02/phasecode/results/enet/sub%02d/fullwindow/sub%02d_enet', subj, subj);
if isempty(tBegin)
    filename = [filename, '.mat'];
    save(filename, 'probability', 'trialinfo', 'toilim', 'dataCW_short', 'dataCCW_short')
else
    filename = [filename, sprintf('_%s', num2str(tBegin)), '.mat'];
    save(filename, 'probability', 'trialinfo', 'toilim')
end
