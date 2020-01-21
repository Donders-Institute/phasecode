function [accuracy] = svm_cichy2(subj, data1, data2, rand1, rand2, groupsize, ngroups, rpt, usecluster, filename, doload)
if ~exist('usecluster', 'var'); usecluster=0; end
doload=1;
if usecluster
        qsubfeval(@svm_cichy, subj, [], [], [], [], [], [], [], [], filename, true,  'timreq', 1800, 'memreq', 15*1024^3, 'batchid', sprintf('cichy_%d_%d', subj, rpt));
%         svm_cichy(subj, [], [], rand1, rand2, [], [], rpt, [], filename, true)
else
    if doload
        load(filename)
        load('/project/3011085.02/phasecode/results/tmpcichy', 'dataCW', 'dataCCW')
        data1=dataCW;
        data2=dataCCW;
    end
    time = -0.1:1/200:1.2;
    for rpt=1:nrpt
%         rng('shuffle')
%             rand1 = randperm(size(data1.trial,1));
%             rand2 = randperm(size(data2.trial,1));
    rand1=truerand(size(data1.trial,1)); [~, rand1] = sort(rand1);
    rand2=truerand(size(data2.trial,1)); [~, rand2] = sort(rand2);
    
    tmpCW = data1;
    tmpCCW = data2;
    tmpCW.trial = zeros(ngroups, size(data1.trial,2));
    tmpCCW.trial = zeros(ngroups, size(data2.trial,2));
    idx=1;
    for k=1:ngroups
        tmpCW.trial(k,:) = mean(data1.trial(rand1(idx:idx+groupsize-1),:,find(time==t)),1);
        tmpCCW.trial(k,:) = mean(data2.trial(rand2(idx:idx+groupsize-1),:,find(time==t)),1);
        idx=idx+groupsize;
    end
tmpCW.time = t;
tmpCCW.time = t;
    
    cfg=[];
    cfg.method = 'crossvalidate';
    cfg.mva=  {dml.svm};
    cfg.statistic = {'accuracy'};
    cfg.type= 'nfold';
    cfg.nfolds = ngroups;
    cfg.design = [ones(size(tmpCW.trial,1),1); 2*ones(size(tmpCCW.trial,1),1)];
        rep.stat(rpt) = ft_timelockstatistics(cfg, tmpCW, tmpCCW);
        accuracy(1,rpt) = rep.stat(rpt).statistic.accuracy;
    end
        filename = sprintf('/project/3011085.02/phasecode/results/tmpcichy_%s_%d', num2str(round(t, 3)), ii);
        save(filename, 'accuracy');
end