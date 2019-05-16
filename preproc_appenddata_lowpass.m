function preproc_appenddata_lowpass(subj)
datainfo;
cfg=[];
path = sprintf('/home/electromag/matves/Data/phasecode/meg_clean/subject%02d/', subj);
cd (path)

if subj==1;
    data{1} = load('session1_200_lowpass.mat');
    data{2} = load('session2_200_lowpass.mat');
    data{3} = load('session3_200_lowpass.mat');
    alldata = ft_appenddata(cfg, data{1}.cleandata, data{2}.cleandata, data{3}.cleandata);
elseif subj==2;
    data{1} = load('session1_200_lowpass.mat');
    data{2} = load('session2_200_lowpass.mat');
    data{3} = load('session3_200_lowpass.mat');
    alldata = ft_appenddata(cfg, data{1}.cleandata, data{2}.cleandata, data{3}.cleandata);
elseif subj==3
    data{1} = load('session1_200_lowpass.mat');
    data{2} = load('session12_200_lowpass.mat');
    data{3} = load('session2_200_lowpass.mat');
    data{4} = load('session3_200_lowpass.mat');
    alldata = ft_appenddata(cfg, data{1}.cleandata, data{2}.cleandata, data{3}.cleandata, data{4}.cleandata);
elseif subj==4
    data{1} = load('session1_200_lowpass.mat');
    data{2} = load('session2_200_lowpass.mat');
    data{3} = load('session22_200_lowpass.mat');
    data{4} = load('session3_200_lowpass.mat');
    alldata = ft_appenddata(cfg, data{1}.cleandata, data{2}.cleandata, data{3}.cleandata, data{4}.cleandata);
elseif subj==5;
    data{1} = load('session1_200_lowpass.mat');
    data{2} = load('session2_200_lowpass.mat');
    data{3} = load('session3_200_lowpass.mat');
    alldata = ft_appenddata(cfg, data{1}.cleandata, data{2}.cleandata, data{3}.cleandata);
elseif subj==6;
    data{1} = load('session1_200_lowpass.mat');
    data{2} = load('session2_200_lowpass.mat');
    data{3} = load('session3_200_lowpass.mat');
    alldata = ft_appenddata(cfg, data{1}.cleandata, data{2}.cleandata, data{3}.cleandata);
elseif subj==7;
    data{1} = load('session1_200_lowpass.mat');
    data{2} = load('session2_200_lowpass.mat');
    data{3} = load('session3_200_lowpass.mat');
    alldata = ft_appenddata(cfg, data{1}.cleandata, data{2}.cleandata, data{3}.cleandata);
elseif subj==8;
    data{1} = load('session1_200_lowpass.mat');
    data{2} = load('session2_200_lowpass.mat');
    data{3} = load('session3_200_lowpass.mat');
    alldata = ft_appenddata(cfg, data{1}.cleandata, data{2}.cleandata, data{3}.cleandata);
elseif subj==9;
    data{1} = load('session1_200_lowpass.mat');
    data{2} = load('session2_200_lowpass.mat');
    data{3} = load('session3_200_lowpass.mat');
    data{4} = load('session32_200_lowpass.mat');
    alldata = ft_appenddata(cfg, data{1}.cleandata, data{2}.cleandata, data{3}.cleandata, data{4}.cleandata);
elseif subj==10;
    data{1} = load('session1_200_lowpass.mat');
    data{2} = load('session2_200_lowpass.mat');
    data{3} = load('session22_200_lowpass.mat');
    data{4} = load('session3_200_lowpass.mat');
    alldata = ft_appenddata(cfg, data{1}.cleandata, data{2}.cleandata, data{3}.cleandata, data{4}.cleandata);
end

%% baseline correct
dt=1/200;
% beamerLag = (2/60);
% lag = dt*round(beamerLag/dt);
lag=35/1000;

cfg=[];
cfg.demean = 'yes';
cfg.baselinewindow = [-0.3+lag -0.1+lag];
alldata = ft_preprocessing(cfg,alldata);

save('alldata200_lowpass_40','alldata','-v7.3');
end

