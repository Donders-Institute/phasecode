% execute the following pipeline for subjects 1-10 to analyze the phasecode
% project.

datainfo;
subj=1;

%% preprocessing
% MEG
for ses=subjects(subj).sessions
  phasecode_preprocessing(subj, ses, 1,1)
end

% Anatomy
phasecode_anatomy

%% analysis
analysis_peakfreq(subj);

analysis_dics(subj);

analysis_lcmv(subj);

for f=[4:1:30 32:2:80]
  analysis_phase(subj, f);
end

% use the virtual channel contralateral to the attended hemifield for
% decoding.
contrasts = {'attended', 'unattended'};
rng('shuffle') % initiate seed with time stamp
for c = 1:numel(contrasts)
  for hemi = 1:2
    for f = [4:1:30 32:2:80]
      % observed data
      rngseed = rand(1)*10^6;
      analysis_decoding(subj,contrasts{c},'hemi', hemi,...
        'do_correcttrials', 1, 'do_avgtrials', 1, 'rngseed', rngseed,...
        'do_prewhiten', 2, 'f', f,'nbins', 18, 'groupsize', 5,'nperm', 20)
      
      % random permutations
      for k=1:10
        rngseed = rand(1)*10^6;
        analysis_decoding(subj,contrasts{c},'hemi', hemi,...
          'do_correcttrials', 1, 'do_avgtrials', 1, 'rngseed', rngseed,...
          'do_prewhiten', 2, 'f', f,'nbins', 18, 'groupsize', 5,...
          'randnr', k,'do_randphasebin',1, 'nrandperm',10, 'nperm', 20);
      end
    end
  end
end


%% Use qsubfeval
subj=5;
rng('shuffle')
for hemi=1:2
  mkdir tmp, cd tmp
for f=4:1:30
rngseed = rand(1)*10^6;
qsubfeval(@analysis_decoding, subj,'attended','hemi', hemi,...
'do_correcttrials', 1, 'do_avgtrials', 1, 'rngseed', rngseed,...
'do_prewhiten', 2, 'f', f,'nbins', 18, 'groupsize', 5,...
'randnr', [],'do_randphasebin',0, 'nrandperm',1, 'nperm', 20, 'memreq', 8*1024^3, 'timreq', 3600);
end
end

% random
subj=1;
rng('shuffle')
for hemi=1:2
  mkdir tmp, cd tmp
for f=4:1:30
for k=1:10
rngseed = rand(1)*10^6;
qsubfeval(@analysis_decoding, subj,'attended','hemi', hemi,...
'do_correcttrials', 1, 'do_avgtrials', 1, 'rngseed', rngseed,...
'do_prewhiten', 2, 'f', f,'nbins', 18, 'groupsize', 5,...
'randnr', k,'do_randphasebin',1, 'nrandperm',10, 'nperm', 20, 'memreq', 8*1024^3, 'timreq', 2*3600);
end
end
end