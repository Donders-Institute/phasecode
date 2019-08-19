function stat = analysis_phasebin_svm(data, trl, timewindow, dorand)
if ~exist('dorand'), dorand=false; end

for l=1:size(trl,2)
    [dat, ntrl] = analysis_preparedata(data, trl{l}, timewindow, 'congruent', dorand);
    
    for hemi=1:size(dat,1)
        cfg=[];
        cfg.method= 'crossvalidate';
        cfg.mva=  {dml.svm};
        cfg.statistic = {'accuracy'};
        cfg.type= 'nfold';
        cfg.nfolds = 5;
        cfg.design = [ones(ntrl(hemi),1); 2*ones(ntrl(hemi),1)];
        stat{l, hemi} = ft_timelockstatistics(cfg, dat{hemi,1}, dat{hemi,2});
    end
end