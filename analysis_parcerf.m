function [tlck, maxchan] = analysis_parcerf(dat, source_parc, atlas)

dat_chan = dat; % keep this in order to find erf window

for cnt = 1:numel(dat)
for k=1:numel(dat{cnt}.trial)
       dat{cnt}.trial{k} = source_parc{1}.F * dat{cnt}.trial{k};
end
   dat{cnt}.label = source_parc{cnt}.label;
end
cfg.appenddim = 'rpt';
data = ft_appenddata(cfg, dat{:});

cfg=[];
cfg.preproc.demean = 'yes';
cfg.preproc.baselinewindow = [-0.1 0];
if 1
cfg.appenddim = 'rpt';
data_chan = ft_appenddata(cfg, dat_chan{:});
tlck = ft_timelockanalysis(cfg, data_chan);

cfgp=[];
cfgp.layout = 'CTF275_helmet.mat';
cfgp.parameter = 'avg';
figure; ft_topoplotER(cfgp, tlck);

x = input('what is the ERF window?')
y = input('what is the ERF peak latency?')
else
    % potentially load the erf peak info
end



cfg.keeptrials = 'yes';
tlck = ft_timelockanalysis(cfg, data);

exclude_label = match_str(atlas.parcellationlabel, {'L_???_01', 'L_MEDIAL.WALL_01', 'R_???_01', 'R_MEDIAL.WALL_01'});
selparc = setdiff(1:numel(atlas.parcellationlabel),exclude_label); % hard coded exclusion of midline and ???

x = nearest(tlck.time, x);
y = nearest(tlck.time, y);
peaksign=sign(mean(mean(tlck.trial(:,:,x(1):x(2)),3),1));

tmp = tlck.trial;
tlck.trial = zeros(size(tmp,1), 374, size(tmp,3));

cnt=0;
for k=selparc
    cnt=cnt+1;
    tlck.trial(:,k,:) = tmp(:,cnt,:) * peaksign(cnt);
end
tlck.label = atlas.parcellationlabel;
tlck.brainordinate = atlas;

% determine maximum index for attend left and right trials
for k=1:2
    idx1 = tlck.trialinfo(:,1)==k;
    if k==1
        idx2 = 188:374;
    else
        idx2 = 1:187;
    end
    [~, maxidx(k)] = max(mean(tlck.trial(idx1, idx2, y)));
end
    maxchan = tlck.label(maxidx);
