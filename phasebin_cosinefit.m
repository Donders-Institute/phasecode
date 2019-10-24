% initialize parameters
freqs = 7:10%[4:1:30 32:2:80];
nperm = 50;
nrandperm = 500;
subj=4;
centerphase = [0 1/3 2/3 1 4/3 5/3]*pi-pi;
resultsdir = '/project/3011085.02/phasecode/results/collapsephasebin/';
contrast = 'attended';
hemis=[1];

for h=hemis
  if (strcmp(contrast, 'attended') || strcmp(contrast, 'unattended'))
    hemi = sprintf('hemi%d_', h);
  else
    hemi = [];
  end
  
  % load the decoding results for every frequency and randomization.
  cnt=1;
  for f=freqs
    for ii = 1:nperm
      tmp{cnt,ii} = load([resultsdir, sprintf('%s/sub%02d_decoding_%sf%d_%d', contrast, subj, hemi, f, ii)]);
      tmp{cnt,ii} = mean(tmp{cnt,ii}.accuracy,2);
    end
    
    % now do the same for random permutations
    for ii = 1:nrandperm
      tmp2{cnt,ii} = load([resultsdir, sprintf('%s/sub%02d_decoding_rand_%sf%d_%d', contrast, subj, hemi, f, ii)]);
      tmp2{cnt,ii} = mean(tmp2{cnt,ii}.accuracy,2);
    end
    cnt=cnt+1;
  end
  
  for k=1:size(tmp,1)
    for ii=1:nperm
      A(ii,:,k) = tmp{k,ii};
    end
    for ii=1:nrandperm
      B(ii,:,k) = tmp2{k,ii};
    end
  end
  % clear tmp tmp2
  
  % shape accuracy per phase bin into fieldtrip structure.
  dat{h}=[];
  dat{h}.dimord = 'rpt_chan_freq';
  dat{h}.freq = freqs;
  dat{h}.trial=permute(A,[2 1 3]);
  for k = 1:nperm
    dat{h}.label{k}=sprintf('accuracy%03d',k);
  end
  
  randdat{h} = dat{h};
  randdat{h}.trial = permute(B, [2,1,3]);
  for k = 1:nrandperm
    randdat{h}.label{k}=sprintf('accuracy%03d',k);
  end
  
  design = centerphase;
  tmpdat=reshape(dat{h}.trial,6,[])';
  tmpcfg=[];
  tmpcfg.cosinefit.statistic='complex';
  s{h}=statfun_cosinefit(tmpcfg,tmpdat,design);
  stat{h} = reshape(s{h}.stat, [nperm, numel(freqs)]);
  r{h} = reshape(s{h}.r, [nperm, numel(freqs)]);
  offset{h} = reshape(s{h}.offset, [nperm, numel(freqs)]);
  
  tmpdat=reshape(randdat{h}.trial,6,[])';
  srand{h}=statfun_cosinefit(tmpcfg,tmpdat,design);
  statrand{h} = reshape(srand{h}.stat, [nrandperm, numel(freqs)]);
  rrand{h} = reshape(srand{h}.r, [nrandperm, numel(freqs)]);
  offsetrand{h} = reshape(srand{h}.offset, [nperm, numel(freqs)]);
end
% tmpdat=reshape(randdat.trial,6,[])';
% for k=1:2600
% tmpdat(k,:) = tmpdat(k, randperm(6));
% end
% srand=statfun_cosinefit(tmpcfg,tmpdat,design);
% statrand = reshape(srand.stat, [numrand, numel(freqs)]);
% rrand = reshape(srand.r, [numrand, numel(freqs)]);
% offsetrand = reshape(srand.offset, [numrand, numel(freqs)]);

cnt=1;
for h=hemis
  avg(cnt,:) = squeeze(mean(mean(dat{h}.trial,2),1))';
  avgr(cnt,:) = squeeze(mean(mean(randdat{h}.trial,2),1))';
  cnt=cnt+1;
end
avg=mean(avg,1);
avgr = mean(avgr,1);
statrand = cat(1,statrand{:});
stat = cat(1,stat{:});


figure; subplot(2,2,1);
plot(freqs, avg);
title('accuracy phasebin')
subplot(2,2,3);
plot(freqs, abs(stat)); hold on;
plot(freqs, abs(mean(stat)), 'k', 'LineWidth',2);
title('cosinefit phasebin')

subplot(2,2,2);
plot(freqs, avgr);
title('accuracy randperm phasebin')
subplot(2,2,4);
plot(freqs, abs(statrand)); hold on;
plot(freqs, abs(mean(statrand)), 'k', 'LineWidth',2);
title('cosinefit randperm phasebin')

%% permutation
num=size(stat,1);
datx=[stat.' statrand.'];
P=blkdiag(ones(num,1)./num,ones(num,1)./num);
obsdiff=abs(datx*P)*[1;-1];
randdiff=zeros(52,100);
for k = 1:1000
  randdiff(:,k) = abs(datx(:,randperm(2*num))*P)*[1;-1];
end

figure;plot(freqs,randdiff,'b')
hold on;plot(freqs,obsdiff,'r','linewidth',2)
title('cosine fit wrt random permutation')

% pvalue?
prb_pos   = zeros(size(obsdiff));

nr=1000;
for k=1:nr
  prb_pos = prb_pos + (obsdiff<randdiff(:,k));
end

prb_pos = prb_pos+1;
nr = nr + 1;
prb_pos = prb_pos/nr;

