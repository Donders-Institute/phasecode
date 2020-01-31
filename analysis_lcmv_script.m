datainfo;
if ~exist('do_virtualchan', 'var'), do_virtualchan = false; end
if ~exist('do_lcmvparc', 'var'), do_lcmvparc = false; end

if do_virtualchan
  [data, tlck, F, L] = analysis_lcmv(subj, '3d');
  % index of maximal gamma power increase
  load([projectdir, sprintf('results/freq/sub%02d_dics_gamma', subj)], 'maxidx', 'Tval')
  [virtualchan, F, L] = analysis_virtualchan(tlck, L, F, maxidx);
  
  save([projectdir, sprintf('results/tlck/sub%02d_virtualchan', subj)], 'virtualchan', 'F', 'L')
end

if do_lcmvparc
  [data, tlck, F, L] = analysis_lcmv(subj, '2d');
  load atlas_subparc374_8k.mat
  for k=1:numel(data)
    data{k}.trial = permute(cat(3, data{k}.trial{:}), [3 1 2]);
    data{k}.time = data{k}.time{1};
    source_parc{k} = lcmvparc(data{k}, tlck{k}, F{k}, L{k}, atlas);
  end
  % align over sessions
  for k=2:numel(source_parc)
    for ch=1:numel(source_parc{1}.label)
      r = corr(source_parc{1}.avg(ch,:)', source_parc{k}.avg(ch,:)');
      if r<0
        source_parc{k}.avg(ch,:) = -source_parc{k}.avg(ch,:);
        source_parc{k}.F(ch,:) = -source_parc{k}.F(ch,:);
      end
    end
  end
  save([projectdir, sprintf('results/tlck/sub%02d_sourceparc.mat', subj)], 'source_parc')
end
