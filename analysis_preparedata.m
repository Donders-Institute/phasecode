function [dat, ntrl] = analysis_preparedata(data, trl, timewindow, contrast, dorand)
if ~exist('dorand','var'); dorand=0; end

cfg=[];
cfg.trials = trl;
% cfg.channel = {'MLO', 'MZO', 'MRO'};
% cfg.preproc.lpfilter = 'yes';
% cfg.preproc.lpfreq = 40;
% cfg.preproc.lpfilttype = 'firws';
trldata = ft_selectdata(cfg, data);
% trldata.trial = smoothdata(trldata.trial, 3, 'movmean', 4);

cfg=[];
cfg.toilim = [timewindow(1) timewindow(end)];
trldata = ft_redefinetrial(cfg, trldata);

switch contrast
    case 'attended'
        idx1 = [11 12; 11 13]; % CW:  first row corresponds to left, second 2 right hemifield
        idx2 = [13 14; 12 14]; % CCW
    case 'congruent'
        idx1 = 11;
        idx2 = 14;
end

for hemi = 1:size(idx1,1) % 1 (left cue), 2 (right cue)
    % for now, take all trials, also incorrect trials, or invalid cued
    % trials
    switch contrast
        case 'attended'
            idx = (trldata.trialinfo(:,1)==hemi);%& alldata.trialinfo(:,4)==hemi %& alldata.trialinfo(:,7)==1
            cfg=[];
            cfg.trials = idx;
            tmpdata = ft_selectdata(cfg, trldata);
        case 'congruent'
%             cfg=[];
%             cfg.trials = (trldata.trialinfo(:,7)==1) & (trldata.trialinfo(:,1)==trldata.trialinfo(:,4));
%             tmpdata = ft_selectdata(cfg, trldata);
            tmpdata = trldata;
    end
    
    cfg=[];
    cfg.trials = find(ismember(tmpdata.trialinfo(:,2), idx1(hemi,:)) | ismember(tmpdata.trialinfo(:,2), idx2(hemi,:)));
    tmpdata = ft_selectdata(cfg, tmpdata);
    if dorand
        rng('shuffle')
        tmpdata.trialinfo = tmpdata.trialinfo(randperm(size(tmpdata.trialinfo,1)),:);
    end
        
    
    % split up conditions
    cfg=[];
    cfg.trials = ismember(tmpdata.trialinfo(:,2), idx1(hemi,:));
    dat{hemi,1} = ft_selectdata(cfg, tmpdata);
    cfg.trials = ismember(tmpdata.trialinfo(:,2), idx2(hemi,:));
    dat{hemi,2} = ft_selectdata(cfg, tmpdata);
    
    ntrl(hemi) = min([numel(dat{hemi,1}.trial), numel(dat{hemi,2}.trial)]);
    for k=1:2
        cfg=[];
        tmpidx = randperm(numel(dat{hemi,k}.trial));
        cfg.trials = tmpidx(1:ntrl(hemi));
        dat{hemi,k} = ft_selectdata(cfg, dat{hemi,k});
    end
    
    % prewhiten
    cfg=[];
    cfg.covariance = 'yes';
    cfg.removemean = 'no';
    for k=1:2
        cov(k) = ft_timelockanalysis(cfg, dat{hemi,k});
    end
    cov = mean(cat(3,cov(:).cov),3);
    cov_inv = cov^-0.5;
    
    for k=1:numel(dat{hemi,1}.trial)
        for l=1:2
            dat{hemi,l}.trial{k} = transpose(transpose(dat{hemi,l}.trial{k})*cov_inv);
        end
    end
    
end
    
