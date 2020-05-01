function analysis_phasic_modulation_batch(subs, hemis, freqs, parcs)

for h=hemis
  for subj=subs
    for w = parcs
      analysis_phasic_modulation(subj, 'freqs', freqs, 'dorand', false, 'doparc', true, 'whichparc', w, 'hemis', h);
    end
  end
end