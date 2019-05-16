function preproc_artifact_rejection(subj,ses)

datainfo;
% subj=5;ses=1;
%load resampled data
% load(sprintf('/home/electromag/matves/Data/phasecode/meg_clean/subject%02d/session%ddownsample.mat', subj, ses), 'datacleanresampled');

% load ecg components
load(sprintf('/home/electromag/matves/Data/phasecode/artifacts/subject%02d/session%dcomp.mat', subj, ses), 'comp');



% Inspect components
cfg = [];
cfg.channel = {comp.label{1:10}}; % components to be plotted
cfg.layout = 'CTF275.lay'; % specify the layout file that should be used for plotting
cfg.compscale = 'local';
ft_databrowser(cfg, comp);

%% Reject components
datainfo;
subs_comp = subjects(subj).session(ses).icacomp;
cfg=[];
cfg.component = subs_comp;
cleandata = ft_rejectcomponent(cfg, comp, datacleanresampled);

% save
save(sprintf('/home/electromag/matves/Data/phasecode/meg_clean/subject%02d/session%d.mat', subj, ses), 'cleandata');

end