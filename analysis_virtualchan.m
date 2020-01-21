function [virtualchan, F, L] = analysis_virtualchan(tlck, L, F, sourceidx)

% make sure there is no ambiguous sign flip ambiguity
c = [1 2; 1 3; 2 3]; c = c(1:numel(tlck),:);
for h=1:2
  tmpF{h} = cat(3,F{:});
  tmpF{h} = squeeze(tmpF{h}(sourceidx(h),:,:));
  if 0
    cnt=1;
    for ii=1:numel(tlck)
      for jj=1:numel(tlck)
        if jj>ii
          tmp(cnt,h) = corr(tmpF{h}(:,ii), tmpF{h}(:,jj));
          cnt=cnt+1;
        end
      end
    end
    x = find(tmp(:,h)<0);
    if ~isempty(x) & numel(x)>1
      x = x(find(c(x(1),:)==c(x(2),:)));
      tmpF{h}(:,x) = -tmpF{h}(:,x);
    end
  else % do it manually by checking the ERF peak
    figure;
    for k=1:numel(tlck)
      v = keepfields(tlck{k}, {'time', 'dimord','trialinfo', 'sampleinfo'});
      v.label{1} = 'tmp';
      for t = 1:size(tlck{k}.trial,3)
        v.trial(:,1,t) = transpose(tmpF{h}(:,k)'*tlck{k}.trial(:,:,t)');
      end
      subplot(1,numel(tlck),k), ft_singleplotER([], v);
    end
    x = input('which virtualchannels should be flipped?')
    tmpF{h}(:,x) = -tmpF{h}(:,x);
  end
end


for k=1:numel(tlck)
  for hemi=1:2
    tmpvirtualchan{hemi} = keepfields(tlck{k}, {'time', 'dimord','trialinfo', 'sampleinfo'});
    tmpvirtualchan{hemi}.label{1} = 'maxgamma';
    tmpvirtualchan{hemi}.trial = zeros(size(tlck{k}.trial,1),1,size(tlck{k}.trial,3));
    for t = 1:size(tlck{k}.trial,3)
      tmpvirtualchan{hemi}.trial(:,1,t) = transpose(tmpF{hemi}(:,k)'*tlck{k}.trial(:,:,t)');
    end
    
    if k==1
      virtualchan{hemi} = tmpvirtualchan{hemi};
    else
      virtualchan{hemi}.trial = [virtualchan{hemi}.trial; tmpvirtualchan{hemi}.trial];
      virtualchan{hemi}.trialinfo = [virtualchan{hemi}.trialinfo; tmpvirtualchan{hemi}.trialinfo];
      virtualchan{hemi}.sampleinfo = [virtualchan{hemi}.sampleinfo; tmpvirtualchan{hemi}.sampleinfo];
    end
  end
end