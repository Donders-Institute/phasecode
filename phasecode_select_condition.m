function [data, idx, idx_ori] = phasecode_select_condition(data, contrast, hemi)

cfg=[];
switch contrast
  case 'congruent'
    idx_ori(1,:) = 11;
    idx_ori(2,:) = 14;
    idx = ismember(data.trialinfo(:,2), idx_ori);
  case 'attended'
    tmpidx{1} = [11 12; 11 13]; % CW:  first row corresponds to left, second to right hemifield
    tmpidx{2} = [13 14; 12 14];
    idx_ori(1,:) = tmpidx{1}(hemi,:);
    idx_ori(2,:) = tmpidx{2}(hemi,:);
    idx = ismember(data.trialinfo(:,2), idx_ori) & data.trialinfo(:,1)==hemi;
  case 'unattended'
    tmpidx{1} = [11 12; 11 13]; % CW:  first row corresponds to left, second to right hemifield
    tmpidx{2} = [13 14; 12 14];
    idx_ori(1,:) = tmpidx{1}(hemi,:);
    idx_ori(2,:) = tmpidx{2}(hemi,:);
    idx = ismember(data.trialinfo(:,2), idx_ori) & data.trialinfo(:,1)~=hemi;
end
cfg.trials = idx;
data = ft_selectdata(cfg, data);