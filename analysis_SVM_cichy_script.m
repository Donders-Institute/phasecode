subj=4;

datainfo;
for ses=1:3
    %     filename = [datadir, sprintf('3016045.07_matves_%03d_%03d', subj, ses), '/cleandata.mat'];
    filename = [projectdir, sprintf('data/sub%02d/sub%02d-meg%02d/sub%02d-meg%02d_cleandata.mat', subj, subj, ses, subj, ses)];
    tmp{ses} = load(filename, 'data');
    tmp{ses} = tmp{ses}.data;
end
alldata = ft_appenddata([], tmp{:});
clear tmp 

cfg=[];
cfg.baseline = 'yes';
cfg.baselinewindow = [-0.1 0];
data = ft_preprocessing(cfg, alldata);

cfg=[];
cfg.latency = [-0.1 1.2];
data = ft_selectdata(cfg, data);
data.trial = permute(cat(3, data.trial{:}), [3 1 2]);
data.time = data.time{1};
data.dimord = 'rpt_chan_time';
%{
% PCA
cfg=[];
cfg.method = 'pca';
comp = ft_componentanalysis(cfg, data);
for k=1:size(comp.trial,2)
v(k) = sum(std(squeeze(comp.trial(:,k,:)),[], 2));
end
sel = find(v/sum(v)>=0.005);
cfg=[];
cfg.component = setdiff(1:size(comp.trial,2), sel);
data = ft_rejectcomponent(cfg, comp, data);
%}

time = data.time;
idx = find(time==0);
% baseline_std = std(data.trial(:,:,1:idx),[],3);
% data.trial = data.trial./repmat(baseline_std, [1 1 length(time)]);

fs = 200%data.fsample;
window = 0.02;
nsample = window*fs;
data.trial = smoothdata(data.trial, 3, 'movmean', nsample);

% select conditions
trl_idxCCW = (data.trialinfo(:,2)==14) & (data.trialinfo(:,7)==1) & (data.trialinfo(:,1)==data.trialinfo(:,4));
trl_idxCW =  (data.trialinfo(:,2)==11) & (data.trialinfo(:,7)==1) & (data.trialinfo(:,1)==data.trialinfo(:,4));
% hemi = 1;%left
% idx1 = [11 12; 11 13]; % CW:  first row corresponds to left, second 2 right hemifield
% idx2 = [13 14; 12 14];
% trl_idxCW = ismember(data.trialinfo(:,2), idx1(hemi,:));
% trl_idxCCW = ismember(data.trialinfo(:,2), idx2(hemi,:));

cor_val_trlCCW = find(trl_idxCCW)'; % all correct and validly cued trial indices (left cued)
cor_val_trlCW = find(trl_idxCW)';

cfg=[];
cfg.trials = cor_val_trlCW;
dataCW = ft_selectdata(cfg, data);
cfg.trials = cor_val_trlCCW;
dataCCW = ft_selectdata(cfg, data);

groupsize = 10;
ntrials = min([size(dataCW.trial,1), size(dataCCW.trial,1)]);
ngroups = floor(ntrials/groupsize);
nrpt = 100;
seq = randseq(10);
save(sprintf('/project/3011085.02/phasecode/results/cichy/sub%02d_seq',subj), 'seq')
if ~exist('usecluster'); usecluster=1; end
filename = sprintf('/project/3011085.02/phasecode/results/cichy/tmpcichy_%s', seq);
    save(filename, 'groupsize', 'nrpt', 'ngroups', 'dataCW', 'dataCCW')
for rpt=1:nrpt
    rand1 = randperm(size(dataCW.trial,1)); 
    rand2 = randperm(size(dataCCW.trial,1));
    svm_cichy(subj, dataCW, dataCCW, rand1, rand2, groupsize, ngroups, rpt, usecluster, filename);
end
   
    
    
    
    


