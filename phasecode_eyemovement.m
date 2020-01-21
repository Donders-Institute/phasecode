datainfo;
subj=4;
cfg=[];
cnt=1;
for ses=subjects(subj).validsessions
  filename = [datadir, sprintf('sub%02d/meg%02d/sub%02d-meg%02d_cleandata.mat', subj, ses, subj, ses)];
  tmp{cnt} = load(filename, 'data');
  tmp{cnt} = removefields(tmp{cnt}.data, 'elec');
  cnt=cnt+1;
end
dat = preprocessing_eyedata(subj, tmp);
data = ft_appenddata([], dat{:});

cfg=[];
cfg.channel = {'visAngleX', 'visAngleY'};
data = ft_selectdata(cfg, data);
data_orig=data;


if ~exist('contrast', 'var'), contrast = 'attended'; end
for hemi = [1 2]
[data, idx, idx_ori] = phasecode_select_condition(data_orig, contrast, hemi);

cfg=[];
cfg.keeptrials = 'yes';
for k=1:size(idx_ori,1)
    cfg.trials = ismember(data.trialinfo(:,2), idx_ori(k,:));
    tl{k} = ft_timelockanalysis(cfg, data);
end

cfgs = [];
cfgs.method = 'montecarlo';
cfgs.statistic = 'indepsamplesT';
cfgs.numrandomization = 1000;
cfgs.correctm = 'cluster';
cfgs.design = [ones(1,size(tl{1}.trial,1)), ones(1,size(tl{2}.trial,1))];
for k=1:2
    neighbours(k).label = tl{1}.label{k};
    neighbours(k).neighblabel = tl{1}.label{mod(k,2)+1};
end
cfgs.neighbours = neighbours;
stat{hemi} = ft_timelockstatistics(cfgs, tl{1}, tl{2});


figure;
for k=1:2
subplot(1,2,k), plot(tl{1}.time, squeeze(mean(tl{1}.trial(:,k,:))));
hold on, plot(tl{2}.time, squeeze(mean(tl{2}.trial(:,k,:))));
end
suptitle('mean eye position for two conditions. X position (left), and Y position (right)')

end