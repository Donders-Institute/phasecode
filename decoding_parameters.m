%% do analyses
subj=4; f=10;
rng('shuffle')
for nbins = [6 12 18]
  for nperm = [10 20]
    mkdir tmp, cd tmp
    for groupsize = [5 10 20]
      filename = '/project/3011085.02/phasecode/results/collapsephasebin/attended/';
      filename = [filename sprintf('%dbins_%dperm_%dgroupsize', nbins, nperm, groupsize)];
      
      for k=1:10
      rngseed = rand(1)*10^6;
        tmpfilename = [filename, sprintf('_%d', k)];
        qsubfeval(@analysis_decoding, subj,'attended','hemi', 2,...
          'do_correcttrials', 1, 'do_avgtrials', 1, 'rngseed',rngseed,...
          'do_prewhiten', 2, 'f', f,'nbins', nbins, 'groupsize', groupsize,...
          'randnr', k,'do_randphasebin',0, 'nrandperm',1, 'nperm', nperm, 'tmpfilename', tmpfilename, 'memreq', 8*1024^3, 'timreq', 3600*2*nbins/6*nperm/10);
        
        rngseed = rand(1)*10^6;
        tmpfilename = [tmpfilename, '_rand'];
        qsubfeval(@analysis_decoding, subj,'attended','hemi', 2,...
          'do_correcttrials', 1, 'do_avgtrials', 1, 'rngseed', rngseed,...
          'do_prewhiten', 2, 'f', f,'nbins', nbins, 'groupsize', groupsize,...
          'randnr', k+10,'do_randphasebin',1, 'nrandperm',10, 'nperm', nperm, 'tmpfilename', tmpfilename, 'memreq', 8*1024^3, 'timreq', 3600*2*nbins/6*nperm/10);
      end
    end
  end
end

%% compare parameters
f=10; subj=4;

x=0;
for nbins = [6 12 18]
  x=x+1;
  y=0;
  for nperm = [10 20]
    y=y+1;
    z=0;
    for groupsize = [5 10 20]
      z=z+1;
      filename = '/project/3011085.02/phasecode/results/collapsephasebin/attended/';
      filename = [filename sprintf('%dbins_%dperm_%dgroupsize', nbins, nperm, groupsize)];
      for k=1:10
        tmpfilename = [filename, sprintf('_%d', k)];
        
        acc_all{x,y,z,k} = load(tmpfilename);
        
        tmpfilename = [tmpfilename, '_rand'];
        accrand_all{x,y,z,k} = load(tmpfilename);
      end
    end
  end
end


for k=1:numel(acc_all)
  acc_all{k} = mean(acc_all{k}.accuracy,3);
  try
    accrand_all{k} = mean(accrand_all{k}.accuracy,4);
  catch
    accrand_all{k} = [];
  end
end

for k=1:numel(acc_all)
  acc_all{k} = mean(acc_all{k});
  try
    accrand_all{k} = squeeze(mean(accrand_all{k},2));
  end
end


c{1} = [0:1/3:5/3]*pi;
c{2} = [0:1/6:11/6]*pi;
c{3} = [0:1/9:17/9]*pi;
[s1, s2, s3, s4] = size(acc_all);
for k=1:s4
  for x=1:s1
    for y=1:s2
      for z=1:s3
        ampl(k,x,y,z) = fit_sine_nobins(acc_all{x,y,z,k}', c{x});
        for kk=1:10
          try
            amplrand(k,kk,x,y,z) = fit_sine_nobins(accrand_all{x,y,z,k}(kk,:)', c{x});
          catch
            amplrand(k,kk,x,y,z) = nan;
          end
        end
      end
    end
  end
end
amplrand = reshape(amplrand, [],s1, s2, s3);

ampl_m = squeeze(mean(ampl,1));
ampl_std = squeeze(std(ampl,[],1));

amplrand_m = squeeze(mean(amplrand,1));
amplrand_std = squeeze(std(amplrand,[],1));

dum = [6 12 18];
figure;
cnt=0;
for x=1:s1
  cnt=cnt+1; subplot(3,4,cnt)
  imagesc([1:2], [1:3], squeeze(ampl_m(x,:,:))', [0.01 0.02]), colorbar;
  title(sprintf('ampl %d bins', dum(x))); xlabel('nperm'), ylabel('groupsize')
  cnt=cnt+1; subplot(3,4,cnt)
  imagesc([1:2], [1:3], squeeze(ampl_std(x,:,:))', [0 5*10^-3]), colorbar;
  title(sprintf('std %d bins', dum(x))); xlabel('nperm'), ylabel('groupsize')
  cnt=cnt+1; subplot(3,4,cnt)
  imagesc([1:2], [1:3], squeeze(amplrand_m(x,:,:))', [0.01 0.02]), colorbar;
  title(sprintf('amplrnd %d bins', dum(x))); xlabel('nperm'), ylabel('groupsize')
  cnt=cnt+1; subplot(3,4,cnt)
  imagesc([1:2], [1:3], squeeze(amplrand_std(x,:,:))', [0 5*10^-3]), colorbar;
  title(sprintf('stdrnd %d bins', dum(x))); xlabel('nperm'), ylabel('groupsize')
end

dum2 = [10 20];
dum3 = [5 10 20];
figure;
cnt=1;
for x=1:s1
%   figure;
%   cnt=1;
  for y=1:s2
    for z=1:s3
      subplot(3,6,cnt),
      boxplot([squeeze(ampl(:,x,y,z)); squeeze(amplrand(:,x,y,z))], [ones(10,1); 2*ones(100,1)])
      title(sprintf('%d bins, %d perm, %d groupsize', dum(x), dum2(y), dum3(z)));
      ylim([0 0.03])
      cnt=cnt+1;
    end
  end
%   suptitle(sprintf('%d bins', dum(x)))
end

for x=1:s1
  for y=1:s2
    for z=1:s3
      for k=1:10
        E(k,x,y,z) = sum(ampl(k,x,y,z)>amplrand(:,x,y,z));
      end
    end
  end
end
Em = squeeze(mean(E));
figure
for x=1:3
  subplot(1,3,x);
  imagesc(1:4,1:3,squeeze(Em(x,:,:)), [90 100]);
  xlabel('groupsize'), ylabel('nperm');
  title(sprintf('%d bins', dum(x)))
end


