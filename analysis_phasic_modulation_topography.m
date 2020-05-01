function stat = analysis_phasic_modulation_topography

datainfo;
load atlas_subparc374_8k.mat
exclude_label = match_str(atlas.parcellationlabel, {'L_???_01', 'L_MEDIAL.WALL_01', 'R_???_01', 'R_MEDIAL.WALL_01'});
selparc = setdiff(1:numel(atlas.parcellationlabel),exclude_label); % hard coded exclusion of midline and ???
for k=1:numel(useparc)
  whichparc{k} = find(contains(atlas.parcellationlabel(selparc), useparc{k}));
end
whichparc = unique(cat(1,whichparc{:}));

hemis=[1 2];
allfreqs = 4:20;
f = [4 9 20; 4 10 13]; % selected frequencies
n = numel(valid_subjects);
for h=hemis
  amp{h} = nan(numel(atlas.parcellationlabel),size(f,2), n);
end


for h=hemis
  freqidx = find(ismember(allfreqs, f(h,:)));
for subj=valid_subjects
  filename = [projectdir, 'results/modulation/'];
  for w=1:numel(selparc)
    
    if ismember(w, whichparc)
      tmp = load([projectdir,'results/modulation/', sprintf('sub%02d_phasicmodulation_decoding_parc_%d', subj, w)], 'amp');
      amp{h}(selparc(w), :,subj) = tmp.amp(h,freqidx);
    else
      tmp = load([projectdir,'results/modulation/', sprintf('sub%02d_phasicmodulation_decoding_hemi%d_parc_%d', subj, h, w)], 'amp');
      amp{h}(selparc(w), :,subj) = tmp.amp(1,:);
    end
  end
end

stat{h}.time = f(h,:);
stat{h}.mean = mean(amp{h},3);
stat{h}.std = std(amp{h},[],3);
stat{h}.dimord = 'chan_time';
stat{h}.label = atlas.parcellationlabel;
stat{h}.brainordinate = atlas;

end

    