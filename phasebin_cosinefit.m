
freqs = [4:1:30 32:2:80];
numrand = 50;
subj=4;
centerphase = [0 1/3 2/3 1 4/3 5/3]*pi-pi;
resultsdir = '/project/3011085.02/phasecode/results/collapsephasebin/';
contrast = 'congruent';
hemi=2;
if (strcmp(contrast, 'attended') || strcmp(contrast, 'unattended'))
  hemi = sprintf('hemi%d_', hemi);
else
  hemi = [];
end
  
  

cnt=1;
for f=freqs
  for ii = 1:numrand
    try
    tmp{cnt,ii} = load([resultsdir, sprintf('%s/sub%02d_decoding_%sf%d_%d', contrast, subj, hemi, f, ii)]);
%     tmp{cnt,ii} = load([resultsdir, sprintf('%s/sub%02d_decoding_%sf%d', contrast, subj, hemi, f)]);
    tmp{cnt,ii} = mean(tmp{cnt,ii}.accuracy,2);
    catch
%       [f ii]
    end
    try
      tmp2{cnt,ii} = load([resultsdir, sprintf('%s/noavg/sub%02d_decoding_rand_%sf%d_%d', contrast, subj, hemi, f, ii)]);
%       tmp2{cnt,ii} = load([resultsdir, sprintf('%s/sub%02d_decoding_rand_%sf%d', contrast, subj, hemi, f)]);
      tmp2{cnt,ii} = mean(tmp2{cnt,ii}.accuracy,2);
    catch
      [f ii]
    end
  end
  cnt=cnt+1;
end

for k=1:size(tmp,1)
  for ii=1:numrand
    A(ii,:,k) = tmp{k,ii};
    B(ii,:,k) = tmp2{k,ii};
  end
end
% clear tmp tmp2

dat=[];
dat.dimord = 'rpt_chan_freq';
dat.freq = freqs;
dat.trial=permute(A,[2 1 3]);
for k = 1:numrand
dat.label{k}=sprintf('accuracy%03d',k);
end

randdat = dat;
randdat.trial = permute(B, [2,1,3]);

design = centerphase;
tmpdat=reshape(dat.trial,6,[])';
tmpcfg=[];
tmpcfg.cosinefit.statistic='complex';
s=statfun_cosinefit(tmpcfg,tmpdat,design);
stat = reshape(s.stat, [numrand, numel(freqs)]);
r = reshape(s.r, [numrand, numel(freqs)]);
offset = reshape(s.offset, [numrand, numel(freqs)]);

tmpdat=reshape(randdat.trial,6,[])';
srand=statfun_cosinefit(tmpcfg,tmpdat,design);
statrand = reshape(srand.stat, [numrand, numel(freqs)]);
rrand = reshape(srand.r, [numrand, numel(freqs)]);
offsetrand = reshape(srand.offset, [numrand, numel(freqs)]);
% tmpdat=reshape(randdat.trial,6,[])';
% for k=1:2600
% tmpdat(k,:) = tmpdat(k, randperm(6));
% end
% srand=statfun_cosinefit(tmpcfg,tmpdat,design);
% statrand = reshape(srand.stat, [numrand, numel(freqs)]);
% rrand = reshape(srand.r, [numrand, numel(freqs)]);
% offsetrand = reshape(srand.offset, [numrand, numel(freqs)]);

figure; subplot(2,2,1);
plot(freqs, squeeze(mean(mean(dat.trial,2),1)));
subplot(2,2,3);
plot(freqs, abs(stat)); hold on;
plot(freqs, abs(mean(stat)), 'k', 'LineWidth',2);

subplot(2,2,2);
plot(freqs, squeeze(mean(mean(randdat.trial,2),1)));
subplot(2,2,4);
plot(freqs, abs(statrand)); hold on;
plot(freqs, abs(mean(statrand)), 'k', 'LineWidth',2);

%% permutation

datx=[stat.' statrand.'];
P=blkdiag(ones(numrand,1)./numrand,ones(numrand,1)./numrand);
obsdiff=abs(datx*P)*[1;-1];
randdiff=zeros(52,100);
for k = 1:1000
randdiff(:,k) = abs(datx(:,randperm(2*numrand))*P)*[1;-1];
end

figure;plot(freqs,randdiff,'b')
hold on;plot(freqs,obsdiff,'r','linewidth',2)

% pvalue?
  prb_pos   = zeros(size(obsdiff));
  
  nr=1000;
  for k=1:nr
    prb_pos = prb_pos + (obsdiff<randdiff(:,k));
  end
  
  prb_pos = prb_pos+1;
  nr = nr + 1;
  prb_pos = prb_pos/nr;
  
