
if ~exist('dorand', 'var'), dorand=false; end
class = 'congruent'; %'attended';
if strcmp(class, 'congruent')
    nhemi = 1;
else
    nhemi = 2;
end
% load data
datainfo;
for k=1:numel(subjects(subj).sessions)
    ses = subjects(subj).sessions(k);
    %     filename = [datadir, sprintf('3016045.07_matves_%03d_%03d', subj, ses), '/cleandata.mat'];
    filename = [projectdir, sprintf('data/sub%02d/sub%02d-meg%02d/sub%02d-meg%02d_cleandata.mat', subj, subj, ses, subj, ses)];
    tmp{k} = load(filename, 'data');
    tmp{k} = tmp{k}.data;
end
data = ft_appenddata([], tmp{:});
data_orig=data;

% divide trials into phase bins
centerphase = [0 0.5 1 1.5]*pi;%[acos(1), acos(0), acos(-1)];
[trl, phase, distance, time] = analysis_alphaphase(data, centerphase);

% for every time window, equalize bins
%{
for k=1:size(trl,1)
    for l=1:size(trl,2)
        tmpn(l) = numel(trl{k,l});
    end
    nsamp = min(tmpn);
    clear tmpn
    for l=1:size(trl,2)
        tmpidx = randperm(numel(trl{k,l}));
        trl{k,l} = trl{k,l}(1:nsamp);
    end
end
%}

%% SVM

accuracy = zeros(numel(time), numel(centerphase),nhemi);
twindow = 1/data_orig.fsample%0.020;

fs = data_orig.fsample;
nsample = fs*twindow;
samprange = (nsample-1)/2;
time_short = time((nsample-1)/2 +1:end-(nsample-1)/2);
%
for k=1:numel(time_short)
    k
    % select data for two conditions
    timewindow = time(k:k+nsample-1);
    
    stat = analysis_phasebin_svm(data_orig, trl(k,:), timewindow, dorand);
    for l=1:numel(centerphase)
        for hemi=size(stat,2)
            accuracy(k,l,hemi) = stat{l, hemi}.statistic.accuracy;
        end
    end
    clear stat;
end
%}

if dorand
    save(sprintf('/project/3011085.02/phasecode/results/phasebin_svm/sub%02d_phasebin_svm1_%d.mat',subj, randnr), 'accuracy');
else
    save(sprintf('/project/3011085.02/phasecode/results/phasebin_svm/sub%02d_phasebin_svm1.mat', subj), 'accuracy');
end

