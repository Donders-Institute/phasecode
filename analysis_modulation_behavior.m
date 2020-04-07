function analysis_modulation_behavior(subj, varargin)
% model whether behavior (reaction times) are modulated by the phase of a
% particular frequency. Use every parcel from a cortical sheet to compute
% the instantaneous phase at stimulus change. Use the phase to bin reaction
% times, and fit a cosine to these bins. Also do this after shuffling
% reaction times (permutation distribution).


contrast = ft_getopt(varargin, 'contrast', 'attended');
method   = ft_getopt(varargin, 'method',   'cosinefit');
model    = ft_getopt(varargin, 'model',    '2d');
freqs    = ft_getopt(varargin, 'freqs',    4:1:30);
nbins    = ft_getopt(varargin, 'nbins',    50);
nperm    = ft_getopt(varargin, 'nperm',    100);
hemis    = ft_getopt(varargin, 'hemis',    [1 2]);
doplot   = ft_getopt(varargin, 'doplot',   false);
dosave   = ft_getopt(varargin, 'dosave',   true);


datainfo;
if strcmp(method, 'kldiv')
  addpath([projectdir, 'scripts/kldiv/']);
end

if strcmp(model, '2d')
  doparc = true;
  load atlas_subparc374_8k.mat
else
  doparc = false;
end
[~, data] = phasecode_getdata(subj, 'doparc', doparc);

% only select validly cued and correct trials
cfg=[];
cfg.trials = (data.trialinfo(:,7)==1) & (data.trialinfo(:,1)==data.trialinfo(:,4));
data = ft_selectdata(cfg, data);

cnt=0;
for f=freqs
  cnt=cnt+1;
  fs = 200;
  % mtmconvol
  numcycles = 2;
  cfg=[];
  cfg.method = 'mtmfft';
  cfg.taper = 'hanning';
  cfg.foi = f;
  cfg.pad = 4;
  cfg.output = 'fourier';
  cfg.t_ftimwin = numcycles*1./cfg.foi;
  cfg2.latency = [1-(numcycles/2*1/f) 1+(numcycles/2*1/f)];
  cfg.keeptrials = 'yes';
  freq = ft_freqanalysis(cfg, ft_selectdata(cfg2, data));
  %     time = cfg.toi;
  %     tcnt(cnt) = nearest(freq.time,1);
  phase{cnt} = angle(squeeze(freq.fourierspctrm))+pi;
end

condition = data.trialinfo(:,1);
rt = data.trialinfo(:,6);
phi=cat(3,phase{:});

bins = [0:2/nbins:2*(nbins-1)/nbins]*pi;
nfreq = numel(freqs);

% make random distributions of RTs
for h=hemis
  trlidx{h} = find(data.trialinfo(:,1)==h);
  tmprt{h} = rt(trlidx{h});
  rtshuff{h} = zeros(numel(tmprt{h}), nperm);
  for k=1:nperm
    rtshuff{h}(:,k) = tmprt{h}(randperm(numel(tmprt{h})));
  end
end

  rtbin = zeros(nfreq, numel(hemis), size(phi,2), nbins);
  rtbinshuff = zeros(nfreq, numel(hemis), size(phi,2), nperm, nbins);
for f=1:nfreq
  f
  for h=hemis
    tmpphi = phi(trlidx{h},:,f);
 
    for k=1:nbins
      centerphase = bins(k);
      lim = [centerphase-pi/2 centerphase+pi/2];
      lim(lim<0) = lim(lim<0) + 2*pi;
      lim(lim>2*pi) = lim(lim>2*pi) - 2*pi;
      if lim(1)>lim(2)
        sel = tmpphi>lim(1) | tmpphi<lim(2);
      else
        sel = tmpphi>lim(1) & tmpphi<lim(2);
      end
      if isempty(sel), break
      end
      for ch = 1:size(sel,2)
        rtbin(f,h,ch,k) = mean(tmprt{h}(sel(:,ch)));
      end
      
      for ii=1:nperm
        for ch = 1:size(sel,2)
          rtbinshuff(f,h,ch,ii,k) = mean(rtshuff{h}(sel(:,ch),ii));
        end
      end
    end
  end
end

switch method
  case 'cosinefit'
    % compute cosine fit
    cfg=[];
    cfg.cosinefit.statistic = 'complex';
    zx = 50;
    centerphase = [0:2/zx:(zx-1)/(zx/2)]*pi-pi;
    
    amp = zeros(nfreq, numel(hemis), size(phi,2));
    ang = zeros(nfreq, numel(hemis), size(phi,2));
    amprand = zeros(nperm, nfreq, numel(hemis), size(phi,2));
    angrand = zeros(nperm, nfreq, numel(hemis), size(phi,2));
    
    [s1, s2, s3, s4] = size(rtbin);
    dat = reshape(rtbin, [], s4);
    s = statfun_cosinefit(cfg, dat, centerphase);
    amp = reshape(abs(s.stat), [s1,s2,s3]);
    ang = reshape(angle(s.stat), [s1,s2,s3]);
    
    [s1, s2, s3, s4, s5] = size(rtbinshuff);
    dat = reshape(rtbinshuff, [], s5);
    s = statfun_cosinefit(cfg, dat, centerphase);
    amprand = reshape(abs(s.stat), [s1,s2,s3,s4]);
    angrand = reshape(angle(s.stat), [s1,s2,s3,s4]);
    
    y = (amp-squeeze(mean(amprand,4)))./squeeze(std(amprand,[],4));
    
  case 'kldiv'
    % compute kl-divergence (deviation from uniformity)
    uniform = ones(nbins,1)/nbins;
    for f=1:nfreq
      for h=hemis
        for ch=1:numel(data.label)
          dat = squeeze(rtbin(f,h,ch,:));
          tmp2 = squeeze(rtbinshuff(f,h,ch,:,:))';
          rtbin_distr = dat/sum(dat);
          kl(f,h,ch) = kldiv([], uniform, rtbin_distr);
          for k=1:nperm
            rtbinshuff_distr = tmp2(:,k)./sum(tmp2(:,k));
            klshuff(k) = kldiv([], uniform, rtbinshuff_distr);
          end
          kl_u(f,h,ch) = mean(klshuff);
          kl_std(f,h,ch) = std(klshuff);
        end
      end
    end
    
    y=(kl-kl_u)./kl_std;
end


stat{1} = rmfield(data, {'trial', 'time'});
stat{1}.dimord ='chan_time';
if strcmp(model, '2d')
  stat{1}.time=freqs;
  exclude_label = match_str(atlas.parcellationlabel, {'L_???_01', 'L_MEDIAL.WALL_01', 'R_???_01', 'R_MEDIAL.WALL_01'});
  selparc = setdiff(1:numel(atlas.parcellationlabel),exclude_label); % hard coded exclusion of midline and ???
  stat{1}.powspctrm = zeros(374,numel(freqs));
  stat{1}.powspctrm(selparc,:) = squeeze(y(:,1,:))';
  stat{1}.brainordinate = atlas;
  load cortex_inflated_shifted.mat
  stat{1}.brainordinate.pos = ctx.pos;
  stat{1}.label = atlas.parcellationlabel;
  stat{2}=stat{1}; stat{2}.powspctrm(selparc,:) = squeeze(y(:,2,:))';
else
  stat{1}.freq=freqs;
  stat{1}.powspctrm = squeeze(y(:,1,:))';
  stat{2}=stat{1}; stat{2}.powspctrm = squeeze(y(:,2,:))';
end

if doplot
  cfgp.layout = 'CTF275_helmet.mat';
  cfgp.colormap = flipud(brewermap(64, 'RdBu'));
  if strcmp(model,'2d')
    cfgp.funparameter = 'powspctrm';
    ft_sourcemovie(cfgp, stat{1})
  else
    figure; ft_topoplotER(cfgp, stat{1}, stat{2})
  end
end

if dosave
  filename = [projectdir, 'results/modulation/', sprintf('sub%02d_cosinefit_behavior_parc.mat', subj)];
  save(filename, 'stat', 'ang','amp', 'amprand', 'angrand','phase', 'rtbin','phi', 'condition')
end
