pp=1;
subs=5:9;
for subj=subs
    
    dt=1/200;
    for jj=1:297
        probabilities{jj} = load(sprintf('/home/electromag/matves/Results/ENET/subject%02d/enet_trial_windowL_%d_%d_MZO.mat', subj, 200, jj),'modelTested', 'ntrials');
    end
    
    
    
    
    probs=[];
    lag=35/1000;
    ntrials=probabilities{1}.ntrials;
    time = (-0.3+dt):dt:(1.2-0.02+dt); %
    % time = time+dt;
    
    for trl = 1:ntrials
        for t=1:length(time)
            probs.trial{trl}(t) = probabilities{t}.modelTested(trl,1);
            probs.trial{ntrials+trl}(t) = probabilities{t}.modelTested(ntrials+trl,2);
        end
        probs.time{trl}=time;
        probs.time{ntrials+trl}=time;
    end
    
    probs.label={'classifier'};
    probs.fsample=200;
    
    %
    
    % probs_tl=ft_timelockanalysis([],probs);
    
    for i=1:ntrials*2;
        alltrials{pp}(i,:)=probs.trial{i}(1,:);
        %     alltrials(i,:)=probs.trial{i}(1,:);
        
    end
    meanprob{pp} = mean(alltrials{pp});
    %%
    
    % cfg=[];
    % cfg.method='mtmfft';
    % cfg.foi=2:1:50;
    % cfg.taper='hanning';
    % fr=ft_freqanalysis(cfg,probs);
    %
    % cfg=[];
    % figure;ft_singleplotER(cfg,fr);
    
    %% EXTRACT
    
%     cfg=[];
%     cfg.latency = [time(1) 0];
%     bl = ft_selectdata(cfg,probs);
    
    % EXTRACT TRIAL
    cfg=[];
    cfg.latency = [0.5+dt 1];
    stim=ft_selectdata(cfg,probs);
    
    %
    cfg=[];
    cfg.method='mtmfft';
    cfg.foi=2:2:50;
    cfg.taper='hanning';
%     cfg.pad = 0.5;
%     bl_fr{pp}=ft_freqanalysis(cfg,bl);
    stim_fr{pp}=ft_freqanalysis(cfg,stim);
    
    cfg=[];
    cfg.parameter = 'powspctrm';
    cfg.operation = 'log10';
    stim_fr_log{pp} = ft_math(cfg,stim_fr{pp});
    % cfg=[];
    % cfg.method='mtmconvol';
    % cfg.foi=2:1:50;
    % cfg.toi=-0.3:0.05:1.2;
    % cfg.taper='hanning';
    % cfg.t_ftimwin=ones(length(cfg.foi))*0.25;
    % cfg.pad = 2;
    % probs_fr_t=ft_freqanalysis(cfg,probs);
    %
    % cfg=[];
    % cfg.baselinetype='relchange';
    % cfg.baseline=[-0.175 -0.125];
    % figure;ft_singleplotTFR(cfg,probs_fr_t);
    %% LOG RATIO
%     cfg=[];
%     cfg.operation='x1/x2';
%     cfg.parameter = 'powspctrm';
%     probability_pow = ft_math(cfg, stim_fr{pp},bl_fr{pp});
%     
%     cfg=[];
%     cfg.operation = 'log10';
%     cfg.parameter = 'powspctrm';
%     logprobability_pow{pp} = ft_math(cfg,probability_pow);
    pp=pp+1;
end

%
% cfg=[];
% cfg.parameter = 'powspctrm';
% GA_logprobability_pow = ft_freqgrandaverage(cfg, logprobability_pow{:});
GA_meanprob = mean(cat(3,meanprob{:}),3);


cfg=[];
cfg.parameter='powspctrm';
ft_singleplotER(cfg,stim_fr_log{:});
avg_stim_fr_log = ft_freqgrandaverage(cfg,stim_fr_log{:});
figure;ft_singleplotER(cfg,avg_stim_fr_log);



%% PLOT
cfg=[];
cfg.parameter = 'powspctrm';
% ft_singleplotER(cfg, stim_pow);
cfg.xlim=[2 50];
figure;
subplot(1,2,2)
ft_singleplotER(cfg,GA_logprobability_pow);
% ft_singleplotER(cfg,logprobability_pow);

xlabel('Frequency (Hz)')
ylabel('Logpower of probability')
title({'powerspectrum stimulus [0.5 1]';'over baseline [-0.3 0]'})
subplot(1,2,1)
plot(time, GA_meanprob)
xlim([time(1) time(end)])
xlabel('Time (s)')
ylabel('Probability')
suptitle({'ENET probability decoding orientation both gratings';sprintf('Subjects 2-9')})


%% ttest?
subj=length(subs);
design = zeros(2,2*subj);
for i = 1:subj
    design(1,i) = i;
end
for i = 1:subj
    design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;

cfg=[];
% cfg.frequency = [];
cfg.method = 'analytic';
cfg.design = design;
cfg.statistic = 'depsamplesT';
cfg.correctm = 'bonferroni';
% cfg.alpha = 0.05;
cfg.parameter = 'powspctrm';
cfg.ivar=2;
cfg.uvar=1;
cfg.frequency = [2 50];
% [stat] = ft_freqstatistics(cfg,stim_fr{:}, stim_fr0{:});
% [stat] = ft_freqstatistics(cfg,logprobability_pow{:}, GA_pow_zero{:});
[stat] = ft_freqstatistics(cfg,stim_fr{:}, bl_fr{:});



cfg=[];
cfg.parameter='stat';
cfg.maskparameter='mask';
figure
ft_singleplotER(cfg,stat);



















%% PHASE ALIGN
clear
ii=1;
classify_on=2;
% hemifield=2;
attention=true;
subs=2:9;


for subj=subs
    % load
    if classify_on==1
        noshift_tmp{ii} = load(sprintf('/home/electromag/matves/Results/SVM/subject%02d/SVM_phase_aligned_orientation_600_noshift.mat', subj));
        shift_tmp{ii} = load(sprintf('/home/electromag/matves/Results/SVM/subject%02d/SVM_phase_aligned_orientation_600.mat', subj));
        ii=ii+1;
    elseif classify_on==2
        for hemifield=1:2
            if attention
                noshift_tmp{ii} = load(sprintf('/home/electromag/matves/Results/SVM/subject%02d/SVM_phase_aligned_hemifield%d_attend_600_noshift.mat', subj, hemifield));
                shift_tmp{ii} = load(sprintf('/home/electromag/matves/Results/SVM/subject%02d/SVM_phase_aligned_hemifield%d_attend_600.mat', subj, hemifield));
            else
                noshift_tmp{ii} = load(sprintf('/home/electromag/matves/Results/SVM/subject%02d/SVM_phase_aligned_hemifield%d_unattend_600_noshift.mat', subj, hemifield));
                shift_tmp{ii} = load(sprintf('/home/electromag/matves/Results/SVM/subject%02d/SVM_phase_aligned_hemifield%d_unattend_600.mat', subj, hemifield));
            end
            ii=ii+1;
        end
    elseif classify_on==3
        noshift_tmp{ii} = load(sprintf('/home/electromag/matves/Results/SVM/subject%02d/SVM_phase_aligned_direction_attention_600_noshift.mat', subj));
        shift_tmp{ii} = load(sprintf('/home/electromag/matves/Results/SVM/subject%02d/SVM_phase_aligned_direction_attention_600.mat', subj));
        ii=ii+1;
    end
end
time = shift_tmp{1}.time;
dt=1/600;

% First make structures with accuracy time series
for subj=1:8
    % shift attend left
    shift{subj}.trial{1} = shift_tmp{subj}.acc';
    shift{subj}.time{1} = time;
    shift{subj}.label = {'classifier'};
    shift{subj}.fsample = 600;
    
    if classify_on==2
        % shift attend right
        shift{subj+8}.trial{1} = shift_tmp{subj+8}.acc';
        shift{subj+8}.time{1} = time;
        shift{subj+8}.label = {'classifier'};
        shift{subj+8}.fsample = 600;
    end
    
    % no shift attend left
    noshift{subj}.trial{1} = noshift_tmp{subj}.acc';
    noshift{subj}.time{1} = time;
    noshift{subj}.label = {'classifier'};
    noshift{subj}.fsample = 600;
    
    if classify_on==2
        % no shift attend right
        noshift{subj+8}.trial{1} = noshift_tmp{subj+8}.acc';
        noshift{subj+8}.time{1} = time;
        noshift{subj+8}.label = {'classifier'};
        noshift{subj+8}.fsample = 600;
    end
end

% now calculate powerspectra
for subj=1:8
    % shift attend left
    cfg=[];
    cfg.latency = [0.5+dt 1];
    window_shift{subj} = ft_selectdata(cfg,shift{subj});
    cfg=[];
    cfg.method='mtmfft';
    cfg.foi=2:2:50;
    cfg.taper='hanning';
    cfg.pad = 0.5;
    shift_pow{subj}=ft_freqanalysis(cfg,window_shift{subj});
    
    if classify_on==2
        % shift attend right
        cfg=[];
        cfg.latency = [0.5+dt 1];
        window_shift{subj+8} = ft_selectdata(cfg,shift{subj+8});
        cfg=[];
        cfg.method='mtmfft';
        cfg.foi=2:2:50;
        cfg.taper='hanning';
        cfg.pad = 0.5;
        shift_pow{subj+8}=ft_freqanalysis(cfg,window_shift{subj+8});
    end
    
    % noshift attend left
    cfg=[];
    cfg.latency = [0.5+dt 1];
    window_noshift{subj} = ft_selectdata(cfg,noshift{subj});
    cfg=[];
    cfg.method='mtmfft';
    cfg.foi=2:2:50;
    cfg.taper='hanning';
    cfg.pad = 0.5;
    noshift_pow{subj}=ft_freqanalysis(cfg,window_noshift{subj});
    
    if classify_on==2
        % noshift attend right
        cfg=[];
        cfg.latency = [0.5+dt 1];
        window_noshift{subj+8} = ft_selectdata(cfg,noshift{subj+8});
        cfg=[];
        cfg.method='mtmfft';
        cfg.foi=2:2:50;
        cfg.taper='hanning';
        cfg.pad = 0.5;
        noshift_pow{subj+8}=ft_freqanalysis(cfg,window_noshift{subj+8});
    end
end

if classify_on==2
    % average powerspectrum over attend left and attend right
    for subj=1:8
        cfg=[];
        avg_shift{subj} = ft_freqgrandaverage(cfg, shift_pow{subj}, shift_pow{subj+8});
        avg_noshift{subj} = ft_freqgrandaverage(cfg, noshift_pow{subj}, noshift_pow{subj+8});
    end
else
    avg_shift = shift_pow;
    avg_noshift = noshift_pow;
end


subj=length(subs);
design = zeros(2,2*subj);
for i = 1:subj
    design(1,i) = i;
end
for i = 1:subj
    design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;

cfg=[];
% cfg.frequency = [];
cfg.method = 'analytic';
cfg.design = design;
cfg.statistic = 'depsamplesT';
cfg.correctm = 'bonferroni';
% cfg.alpha = 0.05;
cfg.parameter = 'powspctrm';
cfg.ivar=2;
cfg.uvar=1;
% [stat] = ft_freqstatistics(cfg,stim_fr{:}, stim_fr0{:});
[stat] = ft_freqstatistics(cfg,avg_shift{:}, avg_noshift{:});

% tmp=stat.mask;
% tmp=double(tmp);
% tmp(tmp==0)=NaN;
% for ii=1:length(tmp)
% significant(2*ii-1) = tmp(ii);
% significant(2*ii) = NaN;
% end
% for ii=1:25
%     if significant(2*ii-1)==1 && ii~=1
%     significant(2*ii-2)=1;
%     significant(2*ii)=1;
%     end
% end

% xlim([stat.time(1) stat.time(end)])
% cfg=[];
% cfg.parameter='stat';
% cfg.maskparameter='mask';
% figure
% ft_singleplotER(cfg,stat);
% title({'decode direction attention'; 't-test SVM logprob realign(10Hz)/normal. vertical lines denote sign value';'powerspectra based on [0.5 1]'})
% xlabel('Frequency (Hz)')
% ylabel('t-value')
% vline(4);vline(34);
%

for subj=1:8
    cfg=[];
    cfg.parameter = 'powspctrm';
    cfg.operation = 'x1/x2';
    pow_shift_noshift{subj} = ft_math(cfg,avg_shift{subj}, avg_noshift{subj});
    cfg=[];
    cfg.parameter ='powspctrm';
    cfg.operation = 'log10';
    logpow_shift_noshift{subj} = ft_math(cfg,pow_shift_noshift{subj});
    avg_shift{subj} = ft_math(cfg,avg_shift{subj});
    avg_noshift{subj} = ft_math(cfg,avg_noshift{subj});
end
GA_shift = ft_freqgrandaverage(cfg, avg_shift{:});
GA_noshift = ft_freqgrandaverage(cfg, avg_noshift{:});
GA_shift_noshift = ft_freqgrandaverage(cfg, logpow_shift_noshift{:});
cfg=[]; cfg.ylim=[-0.6 1];ft_singleplotER(cfg,GA_shift_noshift);
% hold on
% figure(1);plot(significant)
for jj=1:length(subs)
    acc_shift(jj,:) = shift_tmp{jj}.acc';
    acc_noshift(jj,:) = noshift_tmp{jj}.acc';
end
GA_acc_shift = mean(acc_shift);
GA_acc_noshift = mean(acc_noshift);
GA_acc_shift_noshift = GA_acc_shift-GA_acc_noshift;
stat.freq(stat.mask)
%%
% PLOT
figure;
subplot(3,2,2)
cfg=[];
cfg.parameter = 'powspctrm';
ft_singleplotER(cfg,GA_noshift);
xlabel('Frequency (Hz)')
ylabel('Logpower of accuracy')
title('no shift')
subplot(3,2,4)
ft_singleplotER(cfg,GA_shift);
xlabel('Frequency (Hz)')
ylabel('Logpower of accuracy')
title('shift')
subplot(3,2,6)
ft_singleplotER(cfg,GA_shift_noshift);
xlabel('Frequency (Hz)')
ylabel('Logratio of accuracy')
title('shift/noshift')

subplot(3,2,1)
plot(time,smooth(GA_acc_noshift,5))
xlabel('time')
ylabel('classifier accuracy')
title('no shift')
xlim([time(1) time(end)])
subplot(3,2,3)
plot(time,smooth(GA_acc_shift,5))
xlabel('time')
ylabel('classifier accuracy')
title('shift')
xlim([time(1) time(end)])
subplot(3,2,5)
plot(time,smooth(GA_acc_shift_noshift,5))
xlabel('time')
ylabel('difference classifier accuracy')
title('shift - no shift')
xlim([time(1) time(end)])

if classify_on==1
    suptitle({'SVM decoding accuracy with and w/o realignment. decode orientation both gratings';sprintf('Subjects 2-9')})
elseif classify_on==2
    if attention
        suptitle({'SVM decoding accuracy with and w/o realignment. decode orientation attended grating';sprintf('Subjects 2-9')})
    else
        suptitle({'SVM decoding accuracy with and w/o realignment. decode orientation unattended grating';sprintf('Subjects 2-9')})
    end
elseif classify_on==3
    suptitle({'SVM decoding accuracy with and w/o realignment. decode direction of attention';sprintf('Subjects 2-9')})
end






%% treat different conditions as different subjects
clear
ii=1;
for classify_on=1:4
% hemifield=2;
attention=true;
subs=2:9;


for subj=subs
    % load
    if classify_on==1
        noshift_tmp{ii} = load(sprintf('/home/electromag/matves/Results/SVM/subject%02d/SVM_phase_aligned_orientation_600_noshift.mat', subj));
        shift_tmp{ii} = load(sprintf('/home/electromag/matves/Results/SVM/subject%02d/SVM_phase_aligned_orientation_600.mat', subj));
        ii=ii+1;
    elseif classify_on==2
        for hemifield=1:2
                noshift_tmp{ii} = load(sprintf('/home/electromag/matves/Results/SVM/subject%02d/SVM_phase_aligned_hemifield%d_attend_600_noshift.mat', subj, hemifield));
                shift_tmp{ii} = load(sprintf('/home/electromag/matves/Results/SVM/subject%02d/SVM_phase_aligned_hemifield%d_attend_600.mat', subj, hemifield));
            ii=ii+1;
        end
    elseif classify_on==4
        for hemifield=1:2
                noshift_tmp{ii} = load(sprintf('/home/electromag/matves/Results/SVM/subject%02d/SVM_phase_aligned_hemifield%d_unattend_600_noshift.mat', subj, hemifield));
                shift_tmp{ii} = load(sprintf('/home/electromag/matves/Results/SVM/subject%02d/SVM_phase_aligned_hemifield%d_unattend_600.mat', subj, hemifield));
            ii=ii+1;
        end
    elseif classify_on==3
        noshift_tmp{ii} = load(sprintf('/home/electromag/matves/Results/SVM/subject%02d/SVM_phase_aligned_direction_attention_600_noshift.mat', subj));
        shift_tmp{ii} = load(sprintf('/home/electromag/matves/Results/SVM/subject%02d/SVM_phase_aligned_direction_attention_600.mat', subj));
        ii=ii+1;
    end
end
time = shift_tmp{1}.time;
dt=1/600;

% First make structures with accuracy time series
for subj=1:8
    % shift attend left
    shift{subj}.trial{1} = shift_tmp{subj}.acc';
    shift{subj}.time{1} = time;
    shift{subj}.label = {'classifier'};
    shift{subj}.fsample = 600;
    
    if classify_on==2
        % shift attend right
        shift{subj+8}.trial{1} = shift_tmp{subj+8}.acc';
        shift{subj+8}.time{1} = time;
        shift{subj+8}.label = {'classifier'};
        shift{subj+8}.fsample = 600;
    end
    
    % no shift attend left
    noshift{subj}.trial{1} = noshift_tmp{subj}.acc';
    noshift{subj}.time{1} = time;
    noshift{subj}.label = {'classifier'};
    noshift{subj}.fsample = 600;
    
    if classify_on==2
        % no shift attend right
        noshift{subj+8}.trial{1} = noshift_tmp{subj+8}.acc';
        noshift{subj+8}.time{1} = time;
        noshift{subj+8}.label = {'classifier'};
        noshift{subj+8}.fsample = 600;
    end
end

% now calculate powerspectra
for subj=1:8
    % shift attend left
    cfg=[];
    cfg.latency = [0.5+dt 1];
    window_shift{subj} = ft_selectdata(cfg,shift{subj});
    cfg=[];
    cfg.method='mtmfft';
    cfg.foi=2:2:50;
    cfg.taper='hanning';
    cfg.pad = 0.5;
    shift_pow{subj}=ft_freqanalysis(cfg,window_shift{subj});
    
    if classify_on==2
        % shift attend right
        cfg=[];
        cfg.latency = [0.5+dt 1];
        window_shift{subj+8} = ft_selectdata(cfg,shift{subj+8});
        cfg=[];
        cfg.method='mtmfft';
        cfg.foi=2:2:50;
        cfg.taper='hanning';
        cfg.pad = 0.5;
        shift_pow{subj+8}=ft_freqanalysis(cfg,window_shift{subj+8});
    end
    
    % noshift attend left
    cfg=[];
    cfg.latency = [0.5+dt 1];
    window_noshift{subj} = ft_selectdata(cfg,noshift{subj});
    cfg=[];
    cfg.method='mtmfft';
    cfg.foi=2:2:50;
    cfg.taper='hanning';
    cfg.pad = 0.5;
    noshift_pow{subj}=ft_freqanalysis(cfg,window_noshift{subj});
    
    if classify_on==2
        % noshift attend right
        cfg=[];
        cfg.latency = [0.5+dt 1];
        window_noshift{subj+8} = ft_selectdata(cfg,noshift{subj+8});
        cfg=[];
        cfg.method='mtmfft';
        cfg.foi=2:2:50;
        cfg.taper='hanning';
        cfg.pad = 0.5;
        noshift_pow{subj+8}=ft_freqanalysis(cfg,window_noshift{subj+8});
    end
end

if classify_on==1
    jj=0;
elseif classify_on==2
    jj=8;
elseif classify_on==3
    jj=24;
elseif classify_on==4
    jj=16;
end

if classify_on==2
    % average powerspectrum over attend left and attend right
    for subj=1:8
        cfg=[];
        avg_shift{subj+jj} = ft_freqgrandaverage(cfg, shift_pow{subj}, shift_pow{subj+8});
        avg_noshift{subj+jj} = ft_freqgrandaverage(cfg, noshift_pow{subj}, noshift_pow{subj+8});
    end
else
    for subj = 1:8
    avg_shift{subj+jj} = shift_pow{subj};
    avg_noshift{subj+jj} = noshift_pow{subj};
    end
end
end

subj=length(subs);
design = zeros(2,2*subj);
for i = 1:4*subj
    design(1,i) = i;
end
for i = 1:4*subj
    design(1,subj+i) = i;
end
design(2,1:4*subj)        = 1;
design(2,4*subj+1:8*subj) = 2;

cfg=[];
% cfg.frequency = [];
cfg.method = 'analytic';
cfg.design = design;
cfg.statistic = 'depsamplesT';
cfg.correctm = 'no';
% cfg.alpha = 0.05;
cfg.parameter = 'powspctrm';
cfg.ivar=2;
cfg.uvar=1;
% [stat] = ft_freqstatistics(cfg,stim_fr{:}, stim_fr0{:});
[stat] = ft_freqstatistics(cfg,avg_shift{:}, avg_noshift{:});


cfg=[];
cfg.parameter='stat';
cfg.maskparameter='mask';
figure
ft_singleplotER(cfg,stat);
title({'decode direction attention'; 't-test SVM logprob realign(10Hz)/normal. vertical lines denote sign value';'powerspectra based on [0.5 1]'})
xlabel('Frequency (Hz)')
ylabel('t-value')
% vline(4);vline(34);