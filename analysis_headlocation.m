function analysis_headlocation(subj)

if nargin<1 || isempty(subj)
    subj = 4;
end


%% load data, define trials
datainfo;
x = subjects(subj).sessions;
for k=1:numel(x)
cfg=[];
cfg.dataset = subjects(subj).session(x(k)).dataset;
cfg.trialfun = subjects(subj).trialfun;
cfg.continuous = 'yes';
cfg = ft_definetrial(cfg);


%% calculate headposition
% preprocess the headposition data
cfg.channel                 = {'HLC0011','HLC0012','HLC0013', ...
    'HLC0021','HLC0022','HLC0023', ...
    'HLC0031','HLC0032','HLC0033'};
cfg.continuous = 'yes';

headpos = ft_preprocessing(cfg);

% calculate the mean coil position per trial
nTrials = length(headpos.sampleinfo);
for iTrl = 1:nTrials
    coil1(:,iTrl) = [mean(headpos.trial{1,iTrl}(1,:)); mean(headpos.trial{1,iTrl}(2,:)); mean(headpos.trial{1,iTrl}(3,:))];
    coil2(:,iTrl) = [mean(headpos.trial{1,iTrl}(4,:)); mean(headpos.trial{1,iTrl}(5,:)); mean(headpos.trial{1,iTrl}(6,:))];
    coil3(:,iTrl) = [mean(headpos.trial{1,iTrl}(7,:)); mean(headpos.trial{1,iTrl}(8,:)); mean(headpos.trial{1,iTrl}(9,:))];
end

cc{k} = circumcenter(coil1, coil2, coil3);
clear coil*
end

for k=1:numel(cc)
  ntrls(k) = size(cc{k},2);
end
cc = cat(2, cc{:});
cc_mean = mean(cc,2); % mean over trials
cc_dem = [cc - repmat(cc_mean,1,size(cc,2))]';% demean to obtain
% translations and rotations from the average position and orientation
cc_rel = [cc - repmat(cc(:,1),1,size(cc,2))]';% get translations and
% rotations with reference to the first trial
cc=cc';

figure; plot(cc_dem(:,1:3)); 
for k=1:numel(ntrls)-1
  vline(sum(ntrls(1:k)));
end
title('head movement (translation) over sessions'); legend({'x', 'y', 'z'});

figure; plot(cc_dem(:,4:6)); 
for k=1:numel(ntrls)-1
  vline(sum(ntrls(1:k)));
end
title('head rotation over sessions'); legend({'x', 'y', 'z'});

% compute relative change within each session
if numel(x)>3 % MEG system had to be rebooted within 1 session
  ntrls_perses = [0 0 0]';
  for k=1:numel(x)
    tmp = num2str(x(k));
    tmp = str2num(tmp(1));
    ntrls_perses(tmp) = ntrls_perses(tmp) + ntrls(k);
  end
else ntrls_perses = ntrls;
end


idx=1;
for k=1:3
  cctmp = cc(idx:idx+ntrls_perses(k)-1,:);
  ccmean(k,:) = mean(cctmp,1);
  ccstd(k,:) = std(cctmp,[],1);
  idx = idx + ntrls_perses(k);
end
T{1} = 'translation x'; T{2} = 'translation y'; T{3} = 'translation z';
T{4} = 'rotation x'; T{5} = 'rotation y'; T{6} = 'rotation z';
figure
for k=1:6
  subplot(2,3,k);
  errorbar(1:3, ccmean(:,k), ccstd(:,k));
  title(T{k})
  xlim([0 4])
end
suptitle(sprintf('headmotion over sessions, subject %d', subj))
filename = sprintf('/project/3011085.02/phasecode/results/headmotion/sub%02d_headmotion', subj);
saveas(gcf, filename, 'png')
save(fullfile([filename '.mat']), 'cc', 'cc_dem', 'cc_rel', 'ntrls', 'ccmean' ,'ccstd', 'ntrls_perses');
end
