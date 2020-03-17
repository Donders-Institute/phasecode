
datainfo

for subj=valid_subjects
  tmp = load([projectdir, sprintf('results/freq/sub%02d_tfr.mat', subj)], 'tfr_plcmb_high', 'ami_low');
  high{subj} = tmp.tfr_plcmb_high;
  low{subj} = tmp.ami_low;
end

% normalize high frequencies with average over baseline
for k=1:numel(high)
  x1 = nearest(high{k}.time, -0.5);
  x2 = nearest(high{k}.time, -.15);
  high{k}.powspctrm = high{k}.powspctrm./repmat(nanmean(high{k}.powspctrm(:,:,x1:x2),3), [1 1 length(high{k}.time)]) - 1;
end


L = ft_freqgrandaverage([], low{:});
cfg=[];
cfg.appenddim = 'rpt';
L_all = ft_appendfreq(cfg, low{:});
L.std = squeeze(std(L_all.powspctrm));
H = ft_freqgrandaverage([], high{:});
H_all = ft_appendfreq(cfg, high{:});
H.std = squeeze(std(H_all.powspctrm));

save([projectdir, 'results/TFR_group.mat'], 'L', 'H');