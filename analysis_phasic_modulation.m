function analysis_phasic_modulation(subj, varargin)

contrast   = ft_getopt(varargin, 'contrast',   'attended');
hemis      = ft_getopt(varargin, 'hemis',      [1 2]);
freqs      = ft_getopt(varargin, 'freqs',      4:1:30);
npermfiles = ft_getopt(varargin, 'npermfiles', 10);
dorand     = ft_getopt(varargin, 'dorand',     true);
dosave     = ft_getopt(varargin, 'dosave',     true);
doparc     = ft_getopt(varargin, 'doparc',     false);
whichparc  = ft_getopt(varargin, 'whichparc',  []);

datainfo;


resultsdir = [projectdir 'results/collapsephasebin/'];
% loading in data
hcnt=1;
for h=hemis
  if (strcmp(contrast, 'attended') || strcmp(contrast, 'unattended'))
    hemi = sprintf('hemi%d_', h);
  else
    hemi = [];
  end
  
  % load the decoding results for every frequency
  cnt=1;
  for f=freqs
    filedir = [resultsdir, sprintf('%s/sub%02d/', contrast,subj)];
    if doparc
      filedir = [filedir, 'parc/'];
    end
    filename = [filedir, sprintf('sub%02d_decoding_%sf%d', subj, hemi, f)];
    if ~isempty(whichparc)
      filename = [filename, sprintf('_%d', whichparc)];
    end
    tmp{cnt} = load(filename);
    tmp{cnt} = mean(mean(tmp{cnt}.accuracy,3),1);
    
    % now do the same for random phase bins
    if dorand
      for k=1:npermfiles
        if ~isempty(whichparc)
          tmp2{cnt,k} = load([filedir, sprintf('sub%02d_decoding_rand_%sf%d_%d_%d', subj, hemi, f,whichparc, k)]);
        else
          tmp2{cnt,k} = load([filedir, sprintf('sub%02d_decoding_rand_%sf%d_%d', subj, hemi, f, k)]);
        end
        tmp2{cnt,k} = mean(mean(tmp2{cnt,k}.accuracy,4),2);
      end
    end
    cnt=cnt+1;
  end
  if dorand
  if size(tmp2,2)>1
    for k=1:numel(freqs)
      dum = tmp2(k,:);
      tmp2b{k} = squeeze(cat(1,dum{:}));
    end
  else
    for k=1:numel(freqs)
      tmp2b{k} = squeeze(tmp2{k});
    end
  end
  % restructuring decoding accuracies
  acc_rand(hcnt,:,:,:) = cat(3,tmp2b{:});
  end
  acc(hcnt,:,:) = cat(1,tmp{:})';
  
  hcnt=hcnt+1;
end
clear tmp*

zx = size(acc,2);
centerphase = [0:2/zx:(zx-1)/(zx/2)]*pi-pi;


cfg=[];
cfg.cosinefit.statistic = 'complex';
for h=1:hcnt-1
  % observed data
  s = statfun_cosinefit(cfg, squeeze(acc(h,:,:))', centerphase);
  amp(h,:) = abs(s.stat);
  ang(h,:) = angle(s.stat);
  
  % surrogate data
  if dorand
  tmp = reshape(permute(acc_rand, [1 3 4 2]), hcnt-1, length(centerphase), []);
  s = statfun_cosinefit(cfg, squeeze(tmp(h,:,:))', centerphase);
  amp_rand(h,:,:) = reshape(abs(s.stat), numel(freqs),[]);
  ang_rand(h,:,:) = reshape(angle(s.stat), numel(freqs),[]);
  else
    amp_rand=[];
    ang_rand=[];
    acc_rand=[];
  end
end

if dosave
  filedir = [projectdir, 'results/modulation/']; 
  filename = [filedir, sprintf('sub%02d_phasicmodulation_decoding', subj)];
  if numel(hemis)==1
    filename = [filename, sprintf('_hemi%d', hemis)];
  end
  if doparc
    filename = [filename, '_parc'];
    if ~isempty(whichparc)
      filename = [filename, sprintf('_%d', whichparc)];
    end
  end
  save(filename, 'amp', 'ang', 'amp_rand', 'ang_rand', 'acc', 'acc_rand');
end


