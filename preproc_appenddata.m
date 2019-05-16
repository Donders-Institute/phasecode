function preproc_appenddata(subj)
datainfo;
cfg=[];
path = sprintf('/home/electromag/matves/Data/phasecode/meg_clean/subject%02d/', subj);
cd (path)

if subj==1;
    data{1} = load('session1.mat');
    data{2} = load('session2.mat');
    data{3} = load('session3.mat');
    alldata = ft_appenddata(cfg, data{1}.cleandata, data{2}.cleandata, data{3}.cleandata);
elseif subj==2;
    data{1} = load('session1.mat');
    data{2} = load('session2.mat');
    data{3} = load('session3.mat');
    alldata = ft_appenddata(cfg, data{1}.cleandata, data{2}.cleandata, data{3}.cleandata);
elseif subj==3
    data{1} = load('session1.mat');
    data{2} = load('session12.mat');
    data{3} = load('session2.mat');
    data{4} = load('session3.mat');
    alldata = ft_appenddata(cfg, data{1}.cleandata, data{2}.cleandata, data{3}.cleandata, data{4}.cleandata);
elseif subj==4
    data{1} = load('session1.mat');
    data{2} = load('session2.mat');
    data{3} = load('session22.mat');
    data{4} = load('session3.mat');
    alldata = ft_appenddata(cfg, data{1}.cleandata, data{2}.cleandata, data{3}.cleandata, data{4}.cleandata);
elseif subj==5;
    data{1} = load('session1.mat');
    data{2} = load('session2.mat');
    data{3} = load('session3.mat');
    alldata = ft_appenddata(cfg, data{1}.cleandata, data{2}.cleandata, data{3}.cleandata);
elseif subj==6;
    data{1} = load('session1.mat');
    data{2} = load('session2.mat');
    data{3} = load('session3.mat');
    alldata = ft_appenddata(cfg, data{1}.cleandata, data{2}.cleandata, data{3}.cleandata);
elseif subj==7;
    data{1} = load('session1.mat');
    data{2} = load('session2.mat');
    data{3} = load('session3.mat');
    alldata = ft_appenddata(cfg, data{1}.cleandata, data{2}.cleandata, data{3}.cleandata);
elseif subj==8;
    data{1} = load('session1.mat');
    data{2} = load('session2.mat');
    data{3} = load('session3.mat');
    alldata = ft_appenddata(cfg, data{1}.cleandata, data{2}.cleandata, data{3}.cleandata);
elseif subj==9;
    data{1} = load('session1.mat');
    data{2} = load('session2.mat');
    data{3} = load('session3.mat');
    data{4} = load('session32.mat');
    alldata = ft_appenddata(cfg, data{1}.cleandata, data{2}.cleandata, data{3}.cleandata, data{4}.cleandata);
elseif subj==10;
    data{1} = load('session1.mat');
    data{2} = load('session2.mat');
    data{3} = load('session22.mat');
    data{4} = load('session3.mat');
    alldata = ft_appenddata(cfg, data{1}.cleandata, data{2}.cleandata, data{3}.cleandata, data{4}.cleandata);
end


save('alldata','alldata','-v7.3');
end

