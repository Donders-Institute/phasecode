% Test on the group level whether whether behavior (reaction times) is
% modulated by the phase of a particular frequency in a parcel of the
% cortical sheet.
clear 
datainfo;
if ~exist('whichdata', 'var')
    whichdata = input('do you want to test modulation in behavior, or in neural data?');
end
load atlas_subparc374_8k.mat
exclude_label = match_str(atlas.parcellationlabel, {'L_???_01', 'L_MEDIAL.WALL_01', 'R_???_01', 'R_MEDIAL.WALL_01'});
selparc = setdiff(1:numel(atlas.parcellationlabel),exclude_label); % hard coded exclusion of midline and ???

switch whichdata
  case 'behavior'
    freqs = 4:1:30;
  case 'neural'
    freqs = 4:1:20;
end
hemis = [1 2];
nperm = 100;
n = numel(valid_subjects);

for k=1:numel(useparc)
  whichparc{k} = find(contains(atlas.parcellationlabel(selparc), useparc{k}));
end
whichparc = unique(cat(1,whichparc{:}));

for h=hemis
  amp{h} = nan(numel(atlas.parcellationlabel),numel(freqs), n);
  ang{h} = nan(numel(atlas.parcellationlabel),numel(freqs), n);
  ampr{h} = nan(numel(atlas.parcellationlabel),numel(freqs), n, nperm);
end

for subj=1:n
  filename = [projectdir, 'results/modulation/'];
  switch whichdata
    case 'behavior'
      filename = [filename, sprintf('sub%02d_cosinefit_behavior_parc.mat', subj)];
      tmp = load(filename);
      for h=hemis
        amp{h}(selparc, :, subj) = squeeze(tmp.amp(:,h,:))';
        ang{h}(selparc, :, subj) = squeeze(tmp.ang(:,h,:))';
        ampr{h}(selparc, :, subj,:) = permute(squeeze(tmp.amprand(:,h,:,:)),[2 1 3]);
      end
    case 'neural'
      for w=1:numel(whichparc)
        tmp = load([projectdir,'results/modulation/', sprintf('sub%02d_phasicmodulation_decoding_parc_%d', subj, whichparc(w))], 'amp', 'amp_rand', 'ang');
        for h=hemis
          amp{h}(selparc(whichparc(w)), :,subj) = tmp.amp(h,:);
          ampr{h}(selparc(whichparc(w)), :,subj, :) = squeeze(tmp.amp_rand(h,:,:));
          ang{h}(selparc(whichparc(w)), :,subj) = tmp.ang(h,:);
        end
      end
  end
end


% do 2nd level permutation
cfg=[];
cfg.numrandomization = 1000;
cfg.uvar = 2;
cfg.ivar = 1;
cfg.tail = 1;
cfg.clustertail = 1;
cfg.clusteralpha = 0.05;
cfg.alpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.clusterthreshold = 'nonparametric_individual';
cfg.correctm = 'cluster';
cfg.connectivity = eye(374);

if strcmp(whichdata, 'behavior')
tmpamp = amp; % keep the amplitudes outside the ROI for plotting.
tmpampr = ampr;
  % manually select only ROI and FOI
  for h=1:2
    amp{h}(setdiff(1:374, selparc(whichparc)),:,:) = nan;
    amp{h}(:,18:end,:)=[];
    ampr{h}(setdiff(1:374, selparc(whichparc)),:,:,:) = nan;
    ampr{h}(:,18:end,:,:)=[];
  end
  freqs = 4:20;
end

for h=hemis
  [s1 s2 s3 s4] = size(tmpampr{h});
  dat = reshape(amp{h}, s1*s2, s3);
  datrand = reshape(ampr{h}, s1*s2, s3, s4);
  cfg.dim = [s1, s2];
  stat{h} = statfun_perm(cfg, dat, datrand); % s.stat gives mean over subjects. Should this be normalized?
  stat{h}.stat(isnan(stat{h}.stat)) = 0;
  stat{h}.time = freqs;
  stat{h}.dimord = 'chan_time';
  stat{h}.label = atlas.parcellationlabel;
  stat{h}.brainordinate = atlas;
  stat{h}.mean = mean(amp{h},3);
  stat{h}.std = std(amp{h},[],3);
  stat{h}.randstd = mean(std(tmpampr{h}, [], 4),3);
  stat{h}.randmean = mean(mean(tmpampr{h},4),3);
end
for h=1:2
  ampr{h} = nanmean(ampr{h},4);
end

switch whichdata
  case 'behavior'
    filename = [projectdir, 'results/stat_phasicmodulation_behavior'];
    save(filename, 'stat', 'amp', 'ang','ampr')
  case 'neural'
    filename = [projectdir, 'results/stat_phasicmodulation_decoding_parc'];
    save(filename, 'stat', 'amp', 'ang','ampr')
end

n=10;
cfg = [];
cfg.method = 'analytic';
cfg.statistic = 'ft_statfun_bayesfactor';
cfg.ivar = 1;
cfg.uvar = 2;
cfg.design = [ones(1,n), 2*ones(1,n); 1:n, 1:n];
for h=hemis
    dat{h} = [];
    dat{h}.label=atlas.parcellationlabel;
    dat{h}.dimord = 'rpt_chan_freq';
    dat{h}.powspctrm = permute(amp{h},[3 1 2]);
    dat{h}.powspctrm(isnan(dat{h}.powspctrm)) = 0;
    dat{h}.freq = 4:20;
    
    datrand{h} = removefields(dat{h}, 'powspctrm');
    datrand{h}.powspctrm = permute(ampr{h}, [3 1 2]);
    datrand{h}.powspctrm(isnan(datrand{h}.powspctrm)) = 0;
    bf_single(h) = ft_freqstatistics(cfg, dat{h}, datrand{h});
    
    % bf for largest cluster
    tmp = reshape(dat{h}.powspctrm, n, []);
    dat{h}.powspctrm = nanmean(tmp(:, stat{h}.posclusterslabelmat==1),2);
    dat{h}.freq = 0;
    dat{h}.label = {'modulation'};
    tmp = reshape(datrand{h}.powspctrm, n, []);
    datrand{h} = removefields(dat{h}, 'powspctrm');
    datrand{h}.powspctrm = nanmean(tmp(:, stat{h}.posclusterslabelmat==1),2);
    bf(h) = ft_freqstatistics(cfg, dat{h}, datrand{h});
end
save(filename, 'bf','bf_single', '-append')