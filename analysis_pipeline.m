% execute the following pipeline for subjects 1-10 to analyze the phasecode
% project.

datainfo;
do_qsubfeval = true;
subj=input('what is the subject number (1-10)?');
%% preprocessing
% MEG
for ses=subjects(subj).sessions
  phasecode_preprocessing(subj, ses, 1,1)
end

% Anatomy
phasecode_anatomy


%% analysis %
freqs = 4:1:30;
hemis = [1 2];

%%%%%%%%%%%%%%%%%%%%%%
%%% basic analyses %%%
%%%%%%%%%%%%%%%%%%%%%%


% Time-frequency analysis
analysis_timfreq(subj)


% virtual channel based on maximal induced gamma power
analysis_peakfreq(subj);

if do_qsubfeval
  qsubfeval(@analysis_dics, subj, 'memreq', 12*1024^3, 'timreq', 3600);
else
  analysis_dics(subj);
end

if do_qsubfeval
  qsubfeval(@execute_pipeline, 'analysis_lcmv_script', subj, {'do_virtualchan', true}, 'timreq', 3600, 'memreq', 12*1024^3);
else
  execute_pipeline('analysis_lcmv_script', subj, {'do_virtualchan', true});
end

% compute phase at virtual channel
for f=freqs
  if do_qsubfeval
    qsubfeval(@analysis_phase, subj, f, 'timreq', 1800, 'memreq', 6*1024^3);
  else
    analysis_phase(subj, f);
  end
end

%%%%%%%%%%%%%%%%%%%%%%
%%% Basic decoding %%%
%%%%%%%%%%%%%%%%%%%%%%
% show that basic decoding based of stimulus orientation is possible.
% frequency is irrelevant here.
contrasts = {'attended', 'unattended'};
rng('shuffle') % initiate seed with time stamp
for c = 1:numel(contrasts)
  for hemi = hemis % attend left (1) and attend right (2) trials
    rngseed = rand(1)*10^6;
    for do_randphasebin=[0 1]
      if do_qsubfeval
        qsubfeval(@analysis_decoding, subj,contrasts{c},'hemi', hemi,...
          'do_correcttrials', 1, 'do_avgtrials', 1, 'rngseed', rngseed,...
          'do_prewhiten', 2, 'f', 4,'nbins', 1, 'groupsize', 5,'nperm', 100,...
          'do_randphasebin', do_randphasebin, 'memreq', 6*1024^3, 'timreq', 3600);
      else
        analysis_decoding(subj,contrasts{c},'hemi', hemi,...
          'do_correcttrials', 1, 'do_avgtrials', 1, 'rngseed', rngseed,...
          'do_prewhiten', 2, 'f', 4,'nbins', 1, 'groupsize', 5,'nperm', 100,...
          'do_randphasebin', do_randphasebin);
      end
    end
  end
end

% now do group analysis: analysis_decoding_group

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Control analysis: decoding on eye tracker data %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% without binning, frequency is irrelevant
contrasts = {'attended', 'unattended'};
nbins_eye = 1;
rng('shuffle')
rngseed = rand(1)*10^6;
for c = 1:numel(contrasts)
  for do_randphasebin=[0 1]
    for hemi=hemis
      if do_qsubfeval
        qsubfeval(@analysis_decoding, subj,constrasts{c},'hemi', hemi,...
          'do_correcttrials', 1, 'do_avgtrials', 1, 'rngseed', rngseed,...
          'do_prewhiten', 2, 'f', 4,'nbins', nbins_eye, 'groupsize', 5, 'chansel', 'eye',...
          'do_randphasebin',do_randphasebin, 'randnr', [], 'nrandperm',1, 'nperm', 100, 'timreq', 3600, ...
          'memreq', 6*1024^3);
      else
        analysis_decoding(subj,constrasts{c},'hemi', hemi,...
          'do_correcttrials', 1, 'do_avgtrials', 1, 'rngseed', rngseed,...
          'do_prewhiten', 2, 'f', 4,'nbins', nbins_eye, 'groupsize', 5, 'chansel', 'eye',...
          'do_randphasebin',do_randphasebin, 'randnr', [], 'nrandperm',1, 'nperm', 100);
      end
    end
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% phasic modulation in decoding %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use the virtual channel contralateral to the attended hemifield for
% decoding.
contrasts = {'attended'};
rng('shuffle') % initiate seed with time stamp
for c = 1:numel(contrasts)
  for hemi = hemis % attend left (1) and attend right (2) trials
    for f = freqs
      % observed data
      rngseed = rand(1)*10^6;
      if do_qsubfeval
        qsubfeval(@analysis_decoding, subj,contrasts{c},'hemi', hemi,...
          'do_correcttrials', 1, 'do_avgtrials', 1, 'rngseed', rngseed,...
          'do_prewhiten', 2, 'f', f,'nbins', 18, 'groupsize', 5,'nperm', 20,...
          'memreq', 8*1024^3, 'timreq', 3600);
      else
        analysis_decoding(subj,contrasts{c},'hemi', hemi,...
          'do_correcttrials', 1, 'do_avgtrials', 1, 'rngseed', rngseed,...
          'do_prewhiten', 2, 'f', f,'nbins', 18, 'groupsize', 5,'nperm', 20);
      end
      
      % random permutations
      for k=1:10
        rngseed = rand(1)*10^6;
        if do_qsubfeval
          qsubfeval(@analysis_decoding, subj,contrasts{c},'hemi', hemi,...
            'do_correcttrials', 1, 'do_avgtrials', 1, 'rngseed', rngseed,...
            'do_prewhiten', 2, 'f', f,'nbins', 18, 'groupsize', 5,'nperm', 20,...
            'do_randphasebin',1, 'randnr', k, 'nrandperm',10, 'timreq', 2*3600, 'memreq', 10*1024^3);
        else
          analysis_decoding(subj,contrasts{c},'hemi', hemi,...
            'do_correcttrials', 1, 'do_avgtrials', 1, 'rngseed', rngseed,...
            'do_prewhiten', 2, 'f', f,'nbins', 18, 'groupsize', 5,'nperm', 20,...
            'do_randphasebin',1, 'randnr', k, 'nrandperm',10);
        end
      end
    end
  end
end

% model phasic modulation of decoding performance
if do_qsubfeval
  qsubfeval(analysis_phasic_modulation, subj, 'timreq', 3600, 'memreq', 10*1024^3);
else
  analysis_phasic_modulation(subj)
end
% now do group analysis: analysis_phasic_modulation_group



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% phasic modulation of reaction times in 2D space %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute spatial filters for cortical sheet parcels
if do_qsubfeval
  qsubfeval(@execute_pipeline, 'analysis_lcmv_script', subj, {'doparc', true}, 'timreq', 3600, 'memreq', 12*1024^3);
else
  execute_pipeline('analysis_lcmv_script', subj, {'do_virtualchan', true});
end

if do_qsubfeval
  qsubfeval(@analysis_modulation_behavior, subj, 'model', '2d', 'method','cosinefit',...
    'freqs', freqs, 'timreq', 4*3600, 'memreq', 12*1024^3);
else
  analysis_modulation_behavior(subj, 'model', '2d', 'method','cosinefit', 'freqs', freqs);
end
% now do group analysis: analysis_modulation_behavior_group. If anything
% comes out, use that/those parcels in which there is an effect, to
% estimate phase. Then do the following.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% phasic modulation in decoding based FEF/IPS phase %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use parcels of IPS and FEF

load atlas_subparc374_8k.mat
exclude_label = match_str(atlas.parcellationlabel, {'L_???_01', 'L_MEDIAL.WALL_01', 'R_???_01', 'R_MEDIAL.WALL_01'});
selparc = setdiff(1:numel(atlas.parcellationlabel),exclude_label); % hard coded exclusion of midline and ???
for k=1:numel(useparc)
  whichparc{k} = find(contains(atlas.parcellationlabel(selparc), useparc{k}));
end
whichparc = unique(cat(1,whichparc{:}));

%%%%%%%%%%%%
% observed %
%%%%%%%%%%%%
for subj=valid_subjects
  for w=1:numel(whichparc)
    for f=4:20
      for hemi=[1 2]
        qsubfeval(@analysis_decoding, subj,'attended','hemi', hemi,...
          'do_correcttrials', 1, 'do_avgtrials', 1, 'rngseed', [],...
          'do_prewhiten', 2, 'f', f,'nbins', 18, 'groupsize', 5,'doparc', 1, 'whichparc',whichparc(w),...
          'randnr', [],'do_randphasebin',0, 'nrandperm',1, 'nperm', 20, 'memreq', 8*1024^3, 'timreq', 2*3600);
      end
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%
% random permutations %
%%%%%%%%%%%%%%%%%%%%%%%
rng('shuffle')
for k=1:10
  mkdir tmp, cd tmp
  for w=1:numel(whichparc)
    for f=4:20
      rngseed = rand(1)*10^6;
      qsubfeval(@qsub_all, whichparc(w),[1 2], k, f, rngseed, 'timreq', 23*3600, 'memreq', 8*1024^3);
    end
  end
end

% if some jobs crashed:
H=[];
for k=1:10
  for w=1:numel(whichparc)
    for f=4:20
      for subj=valid_subjects
        for hemi=1:2
          if ~exist(sprintf('/project/3011085.06/results/collapsephasebin/attended/sub%02d/parc/sub%02d_decoding_rand_hemi%d_f%d_%d_%d.mat', subj, subj, hemi,f, whichparc(w), k))
            H=[H; subj,hemi, f, whichparc(w), k];
          end
        end
      end
    end
  end
end

cnt=1;
for k=1:ceil(size(H,1)./20)
  try
    qsubfeval(@qsub_tmp, H(cnt:cnt+19,:) , 'timreq', 20*3600, 'memreq', 10*1024^3);
  catch
    qsubfeval(@qsub_tmp, H(cnt:end,:) , 'timreq', 6*3600, 'memreq', 10*1024^3);
  end
  cnt=cnt+20;
end

%%%%%%%%%%%%%%%%%%%%%
% phasic modulation %
%%%%%%%%%%%%%%%%%%%%%
for subj=valid_subjects
  for w=1:numel(whichparc)
    if do_qsubfeval
      qsubfeval(@analysis_phasic_modulation, subj, 'freqs', [4:20], 'dorand', true, 'doparc', true, 'whichparc', whichparc(w), 'hemis', [1 2], 'timreq', 1800,'memreq', 5*1024^3);
    else
      analysis_phasic_modulation(subj, 'freqs', [4:20], 'dorand', true, 'doparc', true, 'whichparc', whichparc(w), 'hemis', [1 2]);
    end
  end
end

%% Group analysis
% TFR on the group level
analysis_timfreq_group

% test whether we can decode the orientation of attended/unattended stimuli
analysis_decoding_group

% test phasic modulation of decoding performance
% FIXME: needs adjustment for statistics
analysis_phasic_modulation_group

whichdata = 'behavior';
analysis_phasic_modulation_parc_group

whichdata = 'neural';
analysis_phasic_modulation_parc_group

%% plotting results

phasecode_plotting

