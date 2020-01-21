function [megdata] = preprocessing_eyedata(subj, megdata, dosave)
% Retrieve eyetracker data from the raw MEG data. From the eyedata those
% trials are selected that are present in the processed MEG data. X- and Y-
% gaze positions on the screen (in pixels) are transformed to X- and Y-
% positions in visual degrees, relative to the central fixation point.
% Eyedata is resampled at 600 Hz and time locked both to stimulus onset and
% stimulus change.

if nargin<1 || isempty(subj), subj = 1; end
if nargin<2 || isempty(megdata), error('please give corresponding MEG data as input'); end
if nargin<3 || isempty(dosave), dosave = false; end

%% Load data and define trials
datainfo; % load subject specific info.
subject = subjects(subj);

x1 = subject.sessions;
if numel(x1)>3 % MEG system had to be rebooted within 1 session
    for k=1:numel(x1)
        tmp = num2str(x1(k));
        x2(k) = str2num(tmp(1));
    end
else x2=x1;
end

ntrl = 0;
for ses=1:3
    s = find(x2==ses);
    
    ii=1;
    for ises = x1(s)
        cfg=[];
        cfg = subjects(subj);
        cfg.dataset = subjects(subj).session(ises).dataset;
        cfg.trialfun = subjects(subj).trialfun;
        cfg = ft_definetrial(cfg);
        cfg.continuous = 'yes';
        cfg.channel = {'UADC005', 'UADC006', 'UADC007'};
        data_tmp{ii} = ft_preprocessing(cfg);
        
        data_tmp{ii}.trialinfo(:, end) = data_tmp{ii}.trialinfo(:,end) + sum(ntrl);
        ntrl = [ntrl, size(data_tmp{ii}.trialinfo,1)]; % manually edit trial number
        
        ii=ii+1;
    end
    
    cfg=[];
    cfg.appenddim = 'rpt';
    data = ft_appenddata(cfg, data_tmp{:});
    clear data_tmp
    
    datatmp = data;
    
    
    %% change X and Y positions into visual angle relative to fixation
    tmp1  =rmfield(data, {'trial', 'label'});
    tmp1.label{1} = 'visAngleX';
    tmp2=tmp1;
    tmp2.label{1} = 'visAngleY';
    [tmp1.trial, tmp2.trial] = transform_eyedata(datatmp);
    
    eyedata{ses} = ft_appenddata([], data, tmp1, tmp2);
    
    % select trials from eye data corresponding to trials in meg data
    if ismember(ses, subject.validsessions)
        % first add the trialnumber of the previous sessions to the current
        % session
        ix = find(ses==subject.validsessions);
        if sum(x2==ses)>1 % more than one data set for a session
            ix = find(ses==subject.validsessions); % the megdata for the current session
            
            % find the first trial of the second dataset in this session
            k=2;
            while megdata{ix}.trialinfo(k,end)>megdata{ix}.trialinfo(k-1,end)
                trlidx = k+1;
                k=k+1;
            end
            
            % hard coded: add trialnumber of first dataset to second dataset
            % also add number of trials from all previous sessions
            megdata{ix}.trialinfo(1:trlidx-1,end) = megdata{ix}.trialinfo(1:trlidx-1,end) + sum(ntrl(1:end-2));
            megdata{ix}.trialinfo(trlidx:end,end) = megdata{ix}.trialinfo(trlidx:end,end) + sum(ntrl(1:end-1));
        else
            megdata{ix}.trialinfo(:,end) = megdata{ix}.trialinfo(:,end) + sum(ntrl(1:end-1));
        end
        
        cfg=[];
        cfg.trials = ismember(eyedata{ses}.trialinfo(:,end), megdata{ix}.trialinfo(:,end));
        eyedata{ses} = ft_selectdata(cfg, eyedata{ses});
    end
end

% select valid sessions
eyedata = eyedata(subject.validsessions);

cfg=[];
cfg.resamplefs = 200;
cfg2=[];
cfg2.appenddim = 'chan';
for k=1:numel(eyedata)
eyedata{k} = ft_resampledata(cfg, eyedata{k});
% for some reason the time axes are not the same. replace manually.
eyedata{k}.time = megdata{k}.time;
if ~isequal(megdata{k}.trialinfo, eyedata{k}.trialinfo)
    error('trialinfo not equal')
end
megdata{k} = ft_appenddata(cfg2, megdata{k}, eyedata{k});
end

if dosave
    save([datadir, sprintf('sub%02d/meg%02d/sub%02d-meg%02d_eyedata.mat', subj, ses, subj, ses)]);
end