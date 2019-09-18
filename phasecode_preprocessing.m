function phasecode_preprocessing(subj, ses, selectartifacts, doica)

if ~exist('selectartifacts', 'var'); selectartifacts = true; end
if ~exist('doica', 'var'); doica = true; end

datainfo; %load subject info

x1 = subjects(subj).sessions;
if numel(x1)>3 % MEG system had to be rebooted within 1 session
  for k=1:numel(x1)
    tmp = num2str(x1(k));
    x2(k) = str2num(tmp(1));
  end
end
s = find(x2==ses);

for ises = x1(s)
  %% define trials
  cfg = subjects(subj);
  cfg.dataset = subjects(subj).session(ises).dataset;
  cfg.trialfun = subjects(subj).trialfun;
  cfg = ft_definetrial(cfg);
  
  % read and preprocess data (filtering, baseline correction etc)
  cfg.channel = subjects(subj).channels;
  cfg.continuous = 'yes';
  cfg.demean = 'yes'; %subtract mean from each channel
  cfg.dftfilter = 'yes';
  cfg.dftfreq = [49 51; 99 101; 149 151];
  cfg.usefftfilt = 'yes';
  cfg.hpfilter = 'yes';
  cfg.hpfreq = 0.1;
  cfg.hpfilttype = 'firws';
  dat{ises} = ft_preprocessing(cfg);
  
  if selectartifacts
    %% Show summary
    %{
    cfg          = [];
    cfg.method   = 'summary';
    cfg.channel = {'MEG'};
    cfg.layout = 'CTF275.lay'; %so timecourses and topographies of individual trials can be looked at
    % cfg.alim     = 1e-12;
    dataRejVis        = ft_rejectvisual(cfg,dat{ises});
    % you can check the rejected trial number by typing
    for i=1:length(dataRejVis.cfg.artfctdef.summary.artifact)
      trlind(i)=find(dat{ises}.sampleinfo(:,1)==dataRejVis.cfg.artfctdef.summary.artifact(i));
    end
    disp(trlind);
    dat{ises}.badtrials=trlind;
    % save data data %could be used to save the bad trialnumbers with the raw
    % data.
    %}
    filename = [projectdir, sprintf('results/artifact/sub%02d_ses%d', subj, ises)];
    load(filename)
    
    
    %% Find jump artifacts
    cfg = [];
    
    % channel selection, cutoff and padding
    cfg.artfctdef.zvalue.channel    = 'MEG';
    cfg.artfctdef.zvalue.cutoff     = 20; % artifact threshold
    cfg.artfctdef.zvalue.trlpadding = 0;
    cfg.artfctdef.zvalue.artpadding = 0;
    cfg.artfctdef.zvalue.fltpadding = 0;
    
    % algorithmic parameters
    cfg.artfctdef.zvalue.cumulative    = 'yes';
    cfg.artfctdef.zvalue.medianfilter  = 'yes';
    cfg.artfctdef.zvalue.medianfiltord = 9;
    cfg.artfctdef.zvalue.absdiff       = 'yes';
    
    % make the process interactive
    cfg.artfctdef.zvalue.interactive = 'yes';
    
    [~, artifact_jump] = ft_artifact_zvalue(cfg, dat{ises});
    
    
    %% Find muscle artifacts
    
    cfg            = [];
    
    % channel selection, cutoff and padding
    cfg.artfctdef.zvalue.channel = 'MEG';
    cfg.artfctdef.zvalue.cutoff      = 4; % artifact threshold
    cfg.artfctdef.zvalue.trlpadding  = 0; % pad at both sides of trial till it is ... seconds.
    cfg.artfctdef.zvalue.fltpadding  = 0;
    cfg.artfctdef.zvalue.artpadding  = 0.2; % take 0.2s padding of an artifact (of the part that is above threshold)
    
    % algorithmic parameters
    cfg.artfctdef.zvalue.bpfilter    = 'yes';
    cfg.artfctdef.zvalue.bpfreq      = [110 140];
    cfg.artfctdef.zvalue.bpfiltord   = 9;
    cfg.artfctdef.zvalue.bpfilttype  = 'but';
    cfg.artfctdef.zvalue.hilbert     = 'yes';
    cfg.artfctdef.zvalue.boxcar      = 0.2;
    
    % make the process interactive
    cfg.artfctdef.zvalue.interactive = 'yes';
    
    [~, artifact_muscle] = ft_artifact_zvalue(cfg, dat{ises});
    %}
    
    %% Find EOG artifacts
    %{
% first define trials only in critical period
cfg            = [];
cfg.dataset = subjects(subj).session(ises).dataset;
cfg.trialfun = 'mytrialfun_critical_period';
cfg = ft_definetrial(cfg);
cfg.datafile = subjects(subj).session(ises).dataset;
cfg.headerfile = subjects(subj).session(ises).dataset;
cfg.continuous = 'yes';

% channel selection, cutoff and padding
cfg.artfctdef.zvalue.channel     = {'UADC005','UADC006'}; % Eye tracker data
% cfg.artfctdef.zvalue.channel     = {'EEG...','EEG...'}; % EOG data.
% Should be one of EEG057 - EEG064
cfg.artfctdef.zvalue.cutoff      = 2;
cfg.artfctdef.zvalue.trlpadding  = 0;
cfg.artfctdef.zvalue.artpadding  = 0.1;
cfg.artfctdef.zvalue.fltpadding  = 0;

% algorithmic parameters
cfg.artfctdef.zvalue.bpfilter   = 'yes';
cfg.artfctdef.zvalue.bpfilttype = 'but';
cfg.artfctdef.zvalue.bpfreq     = [1 15];
cfg.artfctdef.zvalue.bpfiltord  = 4;
cfg.artfctdef.zvalue.hilbert    = 'yes';

% feedback
cfg.artfctdef.zvalue.interactive = 'yes';

[~, artifact_EOG] = ft_artifact_zvalue(cfg);
artifact_visual = artifact_EOG;

    %}
    
    % artfctdef.visual = artifact_EOG;
    load([projectdir, sprintf('results/artifact/sub%02d_ses%d', subj, ises)]);
    artfctdef.muscle = artifact_muscle;
    artfctdef.jump = artifact_jump;
    
    % Save artifacts
    filename = [projectdir, sprintf('results/artifact/sub%02d_ses%d_artifact', subj, ises)];
    save(filename, 'artfctdef')
  else
    % when loading previously defined artifacts:
    filename = [projectdir, sprintf('results/artifact/sub%02d_ses%d_artifact', subj, ises)];
    load(filename)
  end
  tmpartfctdef{ises} = artfctdef;
  % remove trials with high variance
  cfg=[];
  cfg.trials = find(~ismember([1:numel(dat{ises}.trial)], artfctdef.badtrial));
  dat{ises} = ft_selectdata(cfg, dat{ises});
  
  % Reject muslce and jump artifacts from complete data
  cfg = [];
  cfg.artfctdef = removefields(artfctdef, {'visual', 'eye'});
  cfg.artfctdef.reject = 'complete'; % remove complete trials
  cfg.artfctdef.crittoilim = [-1 1]; % only remove trials that have an artifact within first second (grating) or just before the grating onset
  dat{ises} = ft_rejectartifact(cfg, dat{ises});
  
end
rm=[];
for k=1:numel(dat)
  if isempty(dat{k})
    rm = [rm, k];
  end
end
dat(rm)=[];
tmpartfctdef(rm)=[];
data = ft_appenddata([], dat{:});
data.grad = dat{1}.grad;

%% Downsample
cfg=[];
cfg.detrend = 'no';
cfg.demean = 'no';
cfg.resamplefs = 200;
dataresampled = ft_resampledata(cfg,data);

% ICA
filename = [projectdir, sprintf('results/artifact/sub%02d_ses%d_ica', subj, ses)];
if doica
  cfg=[];
  cfg.method          = 'fastica';
  cfg.channel         = 'MEG';
  cfg.fastica.numOfIC = 50;
  comp = ft_componentanalysis(cfg, dataresampled);
  save(filename, 'comp')
else
  load(filename)
end
if ~isfield(comp, 'rejectcomponents')
  cfg = [];
  cfg.channel = {comp.label{1:10}}; % components to be plotted
  cfg.layout = 'CTF275_helmet.mat'; % specify the layout file that should be used for plotting
  cfg.compscale = 'local';
  ft_databrowser(cfg, comp);
  comp.rejectcomponents = input('please enter the to be rejected components')
  
  save(filename, 'comp')
end

cfg           = [];
cfg.component = comp.rejectcomponents;
data         = ft_rejectcomponent(cfg, comp, dataresampled);

% Reject trials with eye artifacts from complete data
for k=1:numel(tmpartfctdef)
cfg = [];
cfg.artfctdef = keepfields(tmpartfctdef{k}, {'visual', 'eye'});
cfg.artfctdef.reject = 'complete'; % remove complete trials
cfg.artfctdef.crittoilim = [-1 1]; % only remove trials that have an artifact within first second (grating) or just before the grating onset
data = ft_rejectartifact(cfg, data);
end

filename = [datadir, sprintf('sub%02d/meg%02d/sub%02d-meg%02d_cleandata.mat', subj,ses, subj, ses)];
save(filename, 'data')

