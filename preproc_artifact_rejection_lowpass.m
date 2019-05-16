function preproc_artifact_rejection_lowpass(subj,ses)

% clear
% subj=10;ses=3;
datainfo;
%load resampled data
load(sprintf('/home/electromag/matves/Data/phasecode/meg_clean/subject%02d/session%ddownsample200_lowpass.mat', subj, ses), 'datacleanresampled');

% load ecg components (SAME FOR BOTH DATASETS)
load(sprintf('/home/electromag/matves/Data/phasecode/artifacts/subject%02d/session%dcomp.mat', subj, ses), 'comp');



% Inspect components
cfg = [];
cfg.channel = {comp.label{1:10}}; % components to be plotted
cfg.layout = 'CTF275.lay'; % specify the layout file that should be used for plotting
cfg.compscale = 'local';
ft_databrowser(cfg, comp);

%% Reject components
% use the same components as for the other (TFA) data. 
datainfo;
subs_comp = subjects(subj).session(ses).icacomp;
cfg=[];
cfg.component = subs_comp;
cleandata = ft_rejectcomponent(cfg, comp, datacleanresampled);
% Data from which the components are rejected is already bandpassed [0.1 40] 


% save
save(sprintf('/home/electromag/matves/Data/phasecode/meg_clean/subject%02d/session%d_200_lowpass.mat', subj, ses), 'cleandata');


%implement this in the classifier script. baseline correct just before
% z-scoring (and then timelock analysis).
% cfg.demean = 'yes'; %subtract mean from each channel
% cfg.baselinewindow = [-0.3+lag -0.1+lag]; % baseline correct data
end