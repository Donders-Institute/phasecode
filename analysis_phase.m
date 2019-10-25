function [phasebin, phase, dist] = analysis_phase(subj, f, centerphase)
% This function loads the original virtualchan and for each trial and frequency
% find the angle over time. It categorizes each trial-freq-timepoint in a
% phasebin and saves this info to disk.
if ~exist('centerphase', 'var') || nargin<3 || isempty(centerphase)
  centerphase = [0 1/3 2/3 1 4/3 5/3]*pi;
end

datainfo;
load([projectdir, sprintf('results/tlck/sub%02d_virtualchan', subj)], 'virtualchan')
addpath('/project/3011085.02/phasecode/scripts/CircStat/');

fs = 200;
% mtmconvol
numcycles = 2;
cfg=[];
cfg.method = 'mtmconvol';
cfg.taper = 'hanning';
cfg.foi = f;
cfg.pad = 4;
cfg.output = 'fourier';
cfg.t_ftimwin = numcycles*1./cfg.foi;
cfg.toi = 1./fs:1./fs:1.2-(numcycles/2*1/f)-1/fs;
cfg.keeptrials = 'yes';
freq = ft_freqanalysis(cfg, virtualchan);
time = cfg.toi;

[s1,s2,s3,s4] = size(freq.fourierspctrm);

phase = angle(squeeze(freq.fourierspctrm))+pi;

for l=1:numel(phase)
  [dist(l), phasebin(l)] = min(abs(circ_dist(phase(l), centerphase)));
end
dist = reshape(dist, size(phase));
phasebin = reshape(phasebin, size(phase));

filename = [projectdir, 'results/phase/', sprintf('sub%02d_phase_%d', subj, f(1))];
save(filename, 'centerphase', 'phase', 'phasebin', 'time');

