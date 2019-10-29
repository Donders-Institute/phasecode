function analysis_peakfreq(subj)

%% load data
datainfo;

cnt=1;
for ses=subjects(subj).validsessions
  filename = [datadir, sprintf('sub%02d/meg%02d/sub%02d-meg%02d_cleandata.mat', subj, ses, subj, ses)];
  tmp{cnt} = load(filename, 'data');
  tmp{cnt} = tmp{cnt}.data;
  cnt=cnt+1;
end
data = ft_appenddata([], tmp{:});
data.grad = tmp{1}.grad;
clear tmp
fs=data.fsample;

cfg=[];
cfg.toilim = [-0.6+1/fs 0];
dataBl = ft_redefinetrial(cfg, data);

cfg.toilim = [0.4+1/fs 1];
dataAct = ft_redefinetrial(cfg, data);

%% power
cfg             = [];
cfg.method      = 'mtmfft';
cfg.output      = 'pow';
cfg.taper       = 'hanning';
cfg.foilim      = [2 100];
cfg.pad         = 1;
cfg.keeptrials  = 'no'; % average baseline over trials
cfg.channel     = {'MLO', 'MRO', 'MZO'};
powBl           = ft_freqanalysis(cfg, dataBl);
powAct          = ft_freqanalysis(cfg, dataAct);

cfg=[];
cfg.parameter = 'powspctrm';
cfg.operation = 'x1./x2-1';
powRatio      = ft_math(cfg, powAct, powBl);

avgPow       = mean(powRatio.powspctrm,1);

gamma.range = [30 90];
range = [find(powRatio.freq==gamma.range(1)) find(powRatio.freq==gamma.range(2))];
[gamma.pow_increase idx]  = max(avgPow(range(1):range(2)));
freqs = powRatio.freq(range(1):range(2));
gamma.peakfreq = freqs(idx);

beta.range = [14 30];
range = [find(powRatio.freq==beta.range(1)) find(powRatio.freq==beta.range(2))];
[beta.pow_increase idx]  = min(avgPow(range(1):range(2)));
freqs = powRatio.freq(range(1):range(2));
beta.peakfreq = freqs(idx);

alpha.range = [8 13];
range = [find(powRatio.freq==alpha.range(1)) find(powRatio.freq==alpha.range(2))];
[alpha.pow_increase idx]  = min(avgPow(range(1):range(2)));
freqs = powRatio.freq(range(1):range(2));
alpha.peakfreq = freqs(idx);

theta.range = [4 7];
range = [find(powRatio.freq==theta.range(1)) find(powRatio.freq==theta.range(2))];
[theta.pow_increase idx]  = min(avgPow(range(1):range(2)));
freqs = powRatio.freq(range(1):range(2));
theta.peakfreq = freqs(idx);

save([projectdir, sprintf('results/freq/sub%02d_peakfreq.mat', subj)], 'theta', 'alpha', 'beta', 'gamma', 'powRatio');