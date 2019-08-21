function [phasebin, phase, dist] = analysis_phase(subj, f, centerphase)
% This function loads the original data and for each trial and frequency
% find the angle over time. It categorizes each trial-freq-timepoint in a
% phasebin and saves this info to disk.
if ~exist('centerphase', 'var') || nargin<3
  centerphase = [0 0.5 1 1.5]*pi;
end

datainfo;
  for k=1:numel(subjects(subj).sessions)
    ses = subjects(subj).sessions(k);
    filename = [datadir, sprintf('3016045.07_matves_%03d_%03d', subj, ses), '/cleandata.mat'];
    filename = [projectdir, sprintf('data/sub%02d/sub%02d-meg%02d/sub%02d-meg%02d_cleandata.mat', subj, subj, ses, subj, ses)];
    tmp{k} = load(filename, 'data');
    tmp{k} = tmp{k}.data;
  end
data = ft_appenddata([], tmp{:});
data.grad = tmp{1}.grad;

% mtmconvol
numcycles = 2;
cfg=[];
cfg.method = 'mtmconvol';
cfg.taper = 'hanning';
cfg.foi = f;
cfg.pad = 4;
cfg.output = 'fourier';
cfg.t_ftimwin = numcycles*1./cfg.foi;
cfg.toi = 1./data.fsample:1./data.fsample:1.2;
cfg.keeptrials = 'yes';
freq = ft_freqanalysis(cfg, data);

[s1,s2,s3,s4] = size(freq.fourierspctrm);
phase = zeros(s1,s4,s3);

for k=1:s1
  for l=1:s3
    t2 = find(isnan(squeeze(freq.fourierspctrm(k,1,l,:))), 1);
    [u, s, v] = svd(squeeze(freq.fourierspctrm(k,:,l,1:min([t2-1, s4]))), 'econ');
    phase(k,:,l) = u(:,1)'*squeeze(freq.fourierspctrm(k,:,l,:));
  end
end
phase = angle(phase)+pi;

tmpcenterphase = [centerphase 2*pi];
phasebin = zeros(s1, s4, s3);
dist = zeros(s1, s4, s3);
for l=1:s3
for k=1:s4
    ptmp = squeeze(phase(:,k,l));
    [tmpdist, idx] = min(transpose(abs(ptmp-tmpcenterphase))); % distance to center phases.
    idx(idx==max(idx)) = min(idx); % smallest and largest angles are the same.
    dist(:,k,l) = tmpdist;
    phasebin(:,k,l) = idx;
end
end

filename = [projectdir, 'results/phase/', sprintf('sub%02d_phase_%d', subj, f(1))];
save(filename, 'centerphase', 'phase', 'phasebin');

