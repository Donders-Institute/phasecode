function [accuracy] = svm_cichy(subj, data1, data2, rand1, rand2, groupsize, ngroups, rpt, usecluster, filename, doload)
if ~exist('usecluster', 'var'); usecluster=0; end

if usecluster
        qsubfeval(@svm_cichy, subj, [], [], rand1, rand2, [], [], rpt, [], filename, true,  'timreq', 1800, 'memreq', 15*1024^3, 'batchid', sprintf('cichy_%d_%d', subj, rpt));
%         svm_cichy(subj, [], [], rand1, rand2, [], [], rpt, [], filename, true)
else
    if doload
        load(filename)
        data1=dataCW;
        data2=dataCCW;
    end
    
    tmpCW = data1;
    tmpCCW = data2;
    tmpCW.trial = zeros(ngroups, size(data1.trial,2), size(data1.trial,3));
    tmpCCW.trial = zeros(ngroups, size(data2.trial,2), size(data2.trial,3));
    idx=1;
    for k=1:ngroups
        tmpCW.trial(k,:,:) = mean(data1.trial(rand1(idx:idx+groupsize-1),:,:),1);
        tmpCCW.trial(k,:,:) = mean(data2.trial(rand2(idx:idx+groupsize-1),:,:),1);
        idx=idx+groupsize;
    end
    
    cfg=[];
    cfg.covariance = 'yes';
    cov1 = ft_timelockanalysis(cfg, tmpCW);
    cov2 = ft_timelockanalysis(cfg, tmpCW);
    cov(:,:,1) = cov1.cov;
    cov(:,:,2) = cov2.cov;
    cov = mean(cov,3);
    cov_inv = cov^-0.5;
    
    for k=1:size(tmpCW.trial,3)
        tmpCW.trial(:,:,k) = tmpCW.trial(:,:,k)*cov_inv;
        tmpCCW.trial(:,:,k) = tmpCCW.trial(:,:,k)*cov_inv;
    end
    
    cfg=[];
    cfg.method = 'crossvalidate';
    cfg.mva=  {dml.svm};
    cfg.statistic = {'accuracy'};
    cfg.type= 'nfold';
    cfg.nfolds = ngroups;
    cfg.design = [ones(size(tmpCW.trial,1),1); 2*ones(size(tmpCCW.trial,1),1)];
    accuracy = zeros(1, numel(data1.time));
    for k=1:numel(data1.time)
        cfg.latency = [data1.time(k) data1.time(k)];
        rep.stat(k) = ft_timelockstatistics(cfg, tmpCW, tmpCCW);
        accuracy(1,k) = rep.stat(k).statistic.accuracy;
    end
    
        filename = [filename, sprintf('_%d', rpt)];
        save(filename, 'accuracy');
end