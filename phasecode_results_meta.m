
%% all trials
datainfo; %load subject info
ntrials=zeros(10,1);

for subj=1:10
for ses=1:3
x1 = subjects(subj).sessions;
if numel(x1)>3 % MEG system had to be rebooted within 1 session
  for k=1:numel(x1)
    tmp = num2str(x1(k));
    x2(k) = str2num(tmp(1));
  end
else x2=x1;
end
s = find(x2==ses);

for ises = x1(s)
  % define trials
  cfg=[];
  cfg = subjects(subj);
  cfg.dataset = subjects(subj).session(ises).dataset;
  cfg.trialfun = subjects(subj).trialfun;
  cfg = ft_definetrial(cfg);
ntrials(subj) = ntrials(subj)+size(cfg.trl,1);
end
end
end


%% trials of interest
for subj=1:10
cfg=[];
cnt=1;
for ses=subjects(subj).validsessions
  filename = [datadir, sprintf('sub%02d/meg%02d/sub%02d-meg%02d_cleandata.mat', subj, ses, subj, ses)];
  tmp{cnt} = load(filename, 'data');
  tmp{cnt} = removefields(tmp{cnt}.data, 'elec');
  cnt=cnt+1;
end
data = ft_appenddata([], tmp{:});
cfg=[];
cfg.trials = (data.trialinfo(:,1)==data.trialinfo(:,4));
data = ft_selectdata(cfg, data);
ntrialsOI(subj) = numel(data.trial);

% behavioral performance
cfg=[];
cfg.trials = (data.trialinfo(:,7)==1);
performance(subj) = sum(cfg.trials)./numel(cfg.trials)*100;
ntrialsAnalysis(subj) = sum(cfg.trials);

clear tmp data
end

save([projectdir, 'results/behavior.mat')