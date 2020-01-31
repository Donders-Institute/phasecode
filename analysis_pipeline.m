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

% test phasic modulation of decoding performance
% FIXME: needs adjustment for statistics
if do_qsubfeval
  qsubfeval(analysis_phasic_modulation, subj, 'timreq', 3600, 'memreq', 10*1024^3);
else
  analysis_phasic_modulation(subj)
end
% now do group analysis: analysis_phasic_modulation_group

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Control analysis: decoding on eye tracker data %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIXME: with our without binning? If with, based on which freq?
f=10;
nbins_eye = 1;
% observed data
rng('shuffle')
rngseed = rand(1)*10^6;
if do_qsubfeval
  qsubfeval(@analysis_decoding, subj,'attended','hemi', hemi,...
    'do_correcttrials', 1, 'do_avgtrials', 1, 'rngseed', rngseed,...
    'do_prewhiten', 2, 'f', f,'nbins', nbins_eye, 'groupsize', 5, 'chansel', 'eye',...
    'do_randphasebin',0, 'randnr', [], 'nrandperm',1, 'nperm', 20, 'timreq', 900, ...
    'memreq', 8*1024^3);
else
  analysis_decoding(subj,'attended','hemi', hemi,...
    'do_correcttrials', 1, 'do_avgtrials', 1, 'rngseed', rngseed,...
    'do_prewhiten', 2, 'f', f,'nbins', nbins_eye, 'groupsize', 5, 'chansel', 'eye',...
    'do_randphasebin',0, 'randnr', [], 'nrandperm',1, 'nperm', 20);
end

% randomized data
for k=1:100
  rngseed = rand(1)*10^6;
  if do_qsubfeval
    qsubfeval(@analysis_decoding, subj,'attended','hemi', hemi,...
      'do_correcttrials', 1, 'do_avgtrials', 1, 'rngseed', rngseed,...
      'do_prewhiten', 2, 'f', f,'nbins', nbins_eye, 'groupsize', 5, 'chansel', 'eye',...
      'do_randphasebin',1, 'randnr', k, 'nrandperm',1, 'nperm', 20,'timreq', 900, ...
      'memreq', 8*1024^3);
  else
    analysis_decoding(subj,'attended','hemi', hemi,...
      'do_correcttrials', 1, 'do_avgtrials', 1, 'rngseed', rngseed,...
      'do_prewhiten', 2, 'f', f,'nbins', nbins_eye, 'groupsize', 5, 'chansel', 'eye',...
      'do_randphasebin',1, 'randnr', k, 'nrandperm',1, 'nperm', 20);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% phasic modulation of reaction times in 2D space %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if do_qsubfeval
  qsubfeval(@analysis_modulation_behavior, subj, 'model', '2d', 'method','cosinefit',...
    'freqs', freqs, 'timreq', 4*1024^3, 'memreq', 12*1024^3);
else
  analysis_modulation_behavior(subj, 'model', '2d', 'method','cosinefit', 'freqs', freqs);
end
% now do group analysis: analysis_modulation_behavior_group. If anything
% comes out, use that/those parcels in which there is an effect, to
% estimate phase. Then do the following.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% phasic modulation in decoding based on parcel phase %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load atlas_subparc374_8k.mat
exclude_label = match_str(atlas.parcellationlabel, {'L_???_01', 'L_MEDIAL.WALL_01', 'R_???_01', 'R_MEDIAL.WALL_01'});
selparc = setdiff(1:numel(atlas.parcellationlabel),exclude_label); % hard coded exclusion of midline and ???
whichparc = match_str(useparc, selparc);
rng('shuffle')
rngseed = rand(1)*10^6;
% observed data
if do_qsubfeval
  qsubfeval(@analysis_decoding, subj,'attended','hemi', hemi,...
    'do_correcttrials', 1, 'do_avgtrials', 1, 'rngseed', rngseed,...
    'do_prewhiten', 2, 'f', f,'nbins', 18, 'groupsize', 5,'doparc', 1, 'whichparc',whichparc,...
    'randnr', [],'do_randphasebin',0, 'nrandperm',1, 'nperm', 20, 'memreq', 8*1024^3, 'timreq', 2*3600);
else
  analysis_decoding, subj,'attended','hemi', hemi,...
    'do_correcttrials', 1, 'do_avgtrials', 1, 'rngseed', rngseed,...
    'do_prewhiten', 2, 'f', f,'nbins', 18, 'groupsize', 5,'doparc', 1, 'whichparc',whichparc,...
    'randnr', [],'do_randphasebin',0, 'nrandperm',1, 'nperm', 20,
end

%random data
for k=1:10
  rngseed = rand(1)*10^6;
  if do_qsubfeval
    qsubfeval(@analysis_decoding, subj,'attended','hemi', hemi,...
      'do_correcttrials', 1, 'do_avgtrials', 1, 'rngseed', rngseed,...
      'do_prewhiten', 2, 'f', f,'nbins', 18, 'groupsize', 5,'doparc', 1, 'whichparc',whichparc,...
      'randnr', k,'do_randphasebin',1, 'nrandperm',1, 'nperm', 20, 'memreq', 8*1024^3, 'timreq', 2*3600);
  else
    analysis_decoding, subj,'attended','hemi', hemi,...
      'do_correcttrials', 1, 'do_avgtrials', 1, 'rngseed', rngseed,...
      'do_prewhiten', 2, 'f', f,'nbins', 18, 'groupsize', 5,'doparc', 1, 'whichparc',whichparc,...
      'randnr', k,'do_randphasebin',1, 'nrandperm',1, 'nperm', 20,
  end
end

if do_qsubfeval
  qsubfeval(analysis_phasic_modulation, subj, 'doparc', true, 'timreq', 3600, 'memreq', 10*1024^3);
else
  analysis_phasic_modulation(subj)
end

%% Group analysis

% test phasic modulation of decoding performance
% FIXME: needs adjustment for statistics
analysis_phasic_modulation_group

analysis_modulation_behavior_group

analysis_phasic_modulation_group('doparc', true)

%% plotting results

phasecode_plotting

