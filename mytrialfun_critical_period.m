function [trl] = mytrialfun_critical_period(cfg)

%NOTE: use this only for artifact selection (only search artifacts in
%critical period.
% this function requires the following fields to be specified
% cfg.dataset
% read the header information and the events from the data
hdr   = ft_read_header(cfg.dataset);
event = ft_read_event(cfg.dataset);
cfg.trialdef.eventtype  = 'UPPT001';
cfg.trialdef.eventtyperesp  = 'UPPT002';
cfg.trialdef.eventvalue = [11 12 13 14];
cfg.trialdef.eventvaluerotation = [21 22 23 24];
cfg.trialdef.eventvalueresp.left = 128; % CCW response
cfg.trialdef.eventvalueresp.right = 8; % CW response


lag = 32.5/1000;
% determine the number of samples before and after the trigger
cfg.trialdef.eventduration = 1; % select only trials with stimulus time of 1s. This is excluding rotation and response!
cfg.trialdef.prestim    = 0.3 - lag; % in seconds (0.8s fixation cross and 0.8s cue)
cfg.trialdef.poststim   = 1+lag; %cfg.trialdef.eventduration + lag; % 
cfg.trialdef.resptime   = 1.4;
prestim  = round(cfg.trialdef.prestim  * hdr.Fs);
poststim = round(cfg.trialdef.poststim * hdr.Fs);

%% Select stimulus and response events
stimEvents = event(strcmp({event.type}, 'UPPT001'));
stimValues = [stimEvents.value];
stimSamples = [stimEvents.sample];
respEvents = event(strcmp({event.type}, 'UPPT002'));
respValues = [respEvents.value];
respSamples = [respEvents.sample];


%% select trials
trl=[];
trlcount = 1; % make sure that the trl matrix only contains valid events

for j=1:length(stimEvents)
  % see whether it is a trigger
  
  %it is a trigger, see whether it is a grating onset event which lasts ~1s
  if j~=1 && j~=length(stimEvents); % j-1 and j+1 cannot be defined.
    if ismember(stimValues(j), cfg.trialdef.eventvalue) && ismember(stimValues(j+1), cfg.trialdef.eventvaluerotation) ...
        && (stimSamples(j+1)-stimSamples(j)) > (hdr.Fs-100) && stimSamples(j+1)-stimSamples(j) < (hdr.Fs + 100)
      begsample     = stimSamples(j) - prestim;
      endsample     = stimSamples(j) + poststim;
      offset        = -prestim;
      cue           = stimValues(j-1); % cue precedes the grating onset trigger
      orientation   = stimValues(j);
      rotationTrig  = stimValues(j+1); % rotation trigger follows the grating onset
      
      % find corresponding reponse
      respIdx = find(respSamples > stimSamples(j), 1); % will give the index of the
      % response which follows the grating onset event
      
      if respSamples(respIdx) - stimSamples(j+1) <= (cfg.trialdef.resptime * hdr.Fs) % if response is within time (1.4sec after onset rotation (=j+1))
        response = respValues(respIdx);
        rt = (respSamples(respIdx) - stimSamples(j+1))/hdr.Fs;
      else % miss
        response = 0;
        rt = 1.41;
      end
      
      % change values of cue and response to 1 (left) and 2 (right)
      cue = cue-1;
      %change response to 1 (left) and 2 (right)
      if response == 128
        response = 1; % left, ccw
      elseif response == 8
        response = 2; % right, cw
      end
      
      % change the rotation in seperate direction and place
      if rotationTrig == 21 || rotationTrig == 23
        rotation = 2; % clockwise
      else
        rotation = 1; % ccw
      end
      
      if rotationTrig == 21 || rotationTrig == 22
        rotateGrating = 1; % left grating rotates
      else
        rotateGrating = 2; % right grating rotates
      end
      
      task = rotation == response; %&& response > 0;
      
      trl(trlcount, :) = [begsample endsample offset cue orientation rotation rotateGrating response rt task trlcount];
      trlcount = trlcount+1;
      
    end
  end
end

end
