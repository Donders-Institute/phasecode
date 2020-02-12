% Test on the group level whether whether behavior (reaction times) is
% modulated by the phase of a particular frequency in a parcel of the
% cortical sheet.

datainfo;
load atlas_subparc374_8k.mat
exclude_label = match_str(atlas.parcellationlabel, {'L_???_01', 'L_MEDIAL.WALL_01', 'R_???_01', 'R_MEDIAL.WALL_01'});
selparc = setdiff(1:numel(atlas.parcellationlabel),exclude_label); % hard coded exclusion of midline and ???

freqs = 4:1:30;
hemis = [1 2];
nperm = 100;
n = numel(valid_subjects);

for h=hemis
  amp{h} = zeros(numel(atlas.parcellationlabel),numel(freqs), n);
  ang{h} = zeros(numel(atlas.parcellationlabel),numel(freqs), n);
  ampr{h} = zeros(numel(atlas.parcellationlabel),numel(freqs), n, nperm);
end

for subj=1:n
  filename = [projectdir, 'results/modulation/', sprintf('sub%02d_cosinefit_behavior_parc.mat', subj)];
  tmp = load(filename);
  for h=hemis
    amp{h}(selparc, :, subj) = squeeze(tmp.amp(:,h,:))';
    
    ang{h}(selparc, :, subj) = squeeze(tmp.ang(:,h,:))';
    
    ampr{h}(selparc, :, subj,:) = permute(squeeze(tmp.amprand(:,h,:,:)),[2 1 3]);
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
cfg.connectivity = full(parcellation2connectivity_midline(atlas));


for h=hemis
  [s1 s2 s3 s4] = size(ampr{h});
  dat = reshape(amp{h}, s1*s2, s3);
  datrand = reshape(ampr{h}, s1*s2, s3, s4);
  cfg.dim = [s1, s2];
  stat{h} = statfun_perm(cfg, dat, datrand); % s.stat gives mean over subjects. Should this be normalized?
  stat{h}.time = freqs;
  stat{h}.dimord = 'chan_time';
  stat{h}.label = atlas.parcellationlabel;
  stat{h}.brainordinate = atlas;
end

filename = [projectdir, 'results/stat_phasicmodulation_behavior'];
save(filename, 'stat')
