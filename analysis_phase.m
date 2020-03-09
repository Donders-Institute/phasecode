function [phasebin, phase, dist, time] = analysis_phase(subj, f, centerphase, numcycles, method, dosave, varargin)
% This function loads the original virtualchan and for each trial and frequency
% find the angle over time. It categorizes each trial-freq-timepoint in a
% phasebin and saves this info to disk.
if ~exist('centerphase', 'var') || nargin<3 || isempty(centerphase)
  centerphase = [0 1/3 2/3 1 4/3 5/3]*pi;
end
if ~exist('numcycles', 'var') || nargin<4 || isempty(numcycles)
  numcycles = 2;
end
if ~exist('method', 'var') || nargin<5 || isempty(method)
  method = 'virtualchan';
end
if ~exist('dosave','var') || nargin<6 || isempty(dosave)
  dosave = true;
end
taper = ft_getopt(varargin, 'taper', 'hanning');

addpath('/project/3011085.02/phasecode/scripts/CircStat/');
datainfo;
switch method
  case 'virtualchan'
    load([projectdir, sprintf('results/tlck/sub%02d_virtualchan', subj)], 'virtualchan');
    dat = virtualchan;
  case 'parc'
    whichparc = ft_getopt(varargin, 'whichparc', 'all');
    load([projectdir, sprintf('results/tlck/sub%02d_sourceparc.mat', subj)])  ;
    [~, dat{1}] = phasecode_getdata(subj, 'doparc', true);
    if ~strcmp(whichparc, 'all')
      cfg=[];
      cfg.channel = dat{1}.label{whichparc};
      dat{1} = ft_selectdata(cfg, dat{1});
    end
end


for h=1:numel(dat)
  fs = 200;
  % mtmconvol
  cfg=[];
  cfg.method = 'mtmconvol';
  if strcmp(taper, 'dpss')
    cfg.taper = 'dpss';
    cfg.tapsmofrq = 4;
  else
    cfg.taper = 'hanning';
  end
  cfg.foi = f;
  cfg.pad = 4;
  cfg.output = 'fourier';
  cfg.t_ftimwin = numcycles*1./cfg.foi;
  cfg.toi = 1./fs:1./fs:1.2-(numcycles/2*1/f)-1/fs;
  cfg.keeptrials = 'yes';
  freq = ft_freqanalysis(cfg, dat{h});
  time = cfg.toi;
  
  [s1,s2,s3,s4] = size(freq.fourierspctrm);
  
  phase{h} = angle(squeeze(freq.fourierspctrm))+pi;
  
  if strcmp(method, 'parc') && strcmp(whichparc, 'all')
    dist=[];
    phasebin=[];
  else
    dist{h} = zeros(size(phase{h}));
    phasebin{h} = zeros(size(phase{h}));
    for l=1:numel(phase{h})
      [dist{h}(l), phasebin{h}(l)] = min(abs(circ_dist(phase{h}(l), centerphase)));
    end
    dist{h} = reshape(dist{h}, size(phase{h}));
    phasebin{h} = reshape(phasebin{h}, size(phase{h}));
  end
end

if dosave
  filename = [projectdir, 'results/phase/', sprintf('sub%02d_phase_%d', subj, f(1))];
  save(filename, 'centerphase', 'phase', 'phasebin', 'time');
end

