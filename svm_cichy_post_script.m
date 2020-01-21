

subj=4;
datadir = '/project/3011085.02/phasecode/results/cichy/';
load([datadir, sprintf('sub%02d_seq', subj)]);

load([datadir, sprintf('tmpcichy_%s',seq)])

for rpt=1:100%nrpt
    tmp = load([datadir, sprintf('tmpcichy_%s_%d',seq, rpt)]);
    accuracy(rpt,:) = tmp.accuracy;
end

figure; plot(-0.1:1/200:1.2, mean(accuracy))
hline(0.5)
time = -0.1:1/200:1.2;
idx = find(time==0.4);
mean(mean(accuracy(:, idx:end)))