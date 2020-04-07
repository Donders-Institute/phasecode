

datainfo;
%% TFR
load([projectdir, 'results/TFR_group.mat'])

% low frequencies
cmap = flipud(brewermap(128, 'RdBu'));
cfgp=[];
cfgp.layout = 'CTF275_helmet.mat';
cfgp.zlim = [-.1 .1];
cfgp.xlim = [-.5 -.15];
cfgp.ylim = [8 13];
cfgp.numcontour = 0;
cfgp.gridscale = 250;
cfgp.marker = 'off';
cfgp.colormap = cmap;
figure; ft_topoplotTFR(cfgp, L);
saveas(gcf, [figures_dir, 'TFR_topo_low.eps'], 'epsc')

cfgp=[];
cfgp.zlim = [-.1 .1];
cfgp.numcontour = 0;
cfgp.gridscale = 250;
cfgp.marker = 'off';
cfgp.colormap = cmap;
cfgp.channel = {'MLO11', 'MLO22', 'MLO23', 'MLO31', 'MLO32', 'MLO33'};
figure; ft_singleplotTFR(cfgp, L);
saveas(gcf, [figures_dir, 'TFR_low_left.svg'])
cfgp.channel = {'MRO11', 'MRO22', 'MRO23', 'MRO31', 'MRO32', 'MRO33'};
figure; ft_singleplotTFR(cfgp, L);
saveas(gcf, [figures_dir, 'TFR_low_right.svg'])


% high frequencies

cmap = flipud(brewermap(128, 'RdBu'));
cmap(1:64,:) = [];

chans = {'MLO11', 'MLO21', 'MLO22', 'MLO31', 'MRO11', ...
  'MRO21', 'MRO22', 'MRO31', 'MZO01', 'MZO02'};

cfgp=[];
cfgp.zlim = [0 0.25];
cfgp.colormap = cmap;
cfgp.xlim = [-0.4 1];
cfgp.channel = chans;
figure; ft_singleplotTFR(cfgp, H);
saveas(gcf, [figures_dir, 'TFR_high.svg'])


cfgp=removefields(cfgp, 'channel');
cfgp.layout = 'CTF275_helmet.mat';
cfgp.zlim = [0 0.25];
cfgp.xlim = [0.4 1];
cfgp.ylim = [50 65];
cfgp.numcontour = 0;
cfgp.gridscale = 250;
cfgp.marker = 'off';
cfgp.highlight = 'on';
cfgp.highlightchannel = chans;
cfgp.highlightsymbol = 'O';
cfgp.highlightsize = 8;
cfgp.highlightcolor = [1 1 1];
figure; ft_topoplotTFR(cfgp, H);
saveas(gcf, [figures_dir, 'TFR_topo_high.eps'], 'epsc')

cfg=[];
cfg.avgoverfreq = 'yes';
cfg.avgovertime = 'yes';
cfg.avgoverchan = 'yes';
cfg.latency = [0.2 1];
cfg.frequency = [45 65];
cfg.channel = chans;
H_avg = ft_selectdata(cfg, H);

sprintf('mean gamma power increase %d percent, STD %d percent', H_avg.powspctrm*100, H_avg.std*100)

%% source location gamma virtual channel
% get average power increase
increase=[];
for subj=1:10
  tmp = load([projectdir, sprintf('results/freq/sub%02d_dics_gamma', subj)], 'maxidx', 'powRatio');
  pow(subj,:) = tmp.powRatio.mean-1;
  increase = [increase; 100*pow(subj,tmp.maxidx)];
end
load standard_sourcemodel3d6mm.mat
P = sourcemodel;
P.avg = 100*nanmean(pow)';

% interpolate source onto cortical sheet
m1 = load('surface_white_both.mat');
m2 = load('surface_white_both_shifted');
cfg=[];
cfg.parameter = {'avg'};
sint=ft_sourceinterpolate(cfg, P, m1.mesh);
sint.hemisphere = m2.mesh.hemisphere;
loc=[];
for subj=1:10
tmp = load([projectdir, sprintf('results/freq/sub%02d_dics_gamma', subj)], 'Tval');
tmp.Tval.pos = sourcemodel.pos;
loc(subj,:) = find_interpolated_location(tmp.Tval, 'stat', sint);
end

% plot
cfgp=[];
cfgp.method = 'surface';
cmap = flipud(brewermap(128, 'RdBu'));
cmap(1:64,:) = [];
cfgp.funcolormap = cmap;
cfgp.funparameter = 'avg';
cfgp.funcolorlim = [0 80];
cfgp.colorbar = 'no';
cfgp.surfinflated = 'surface_white_both_shifted.mat';
ft_sourceplot(cfgp, sint)

virtualchanpos=[];
virtualchanpos.chanpos = m2.mesh.pos(loc,:);
virtualchanpos.elecpos = m2.mesh.pos(loc,:);
for k=1:numel(loc)
virtualchanpos.label{k}   = sprintf('%d', k);
end
ft_plot_sens(virtualchanpos,'elecshape','sphere', 'elecsize', 8, 'facecolor', 'black')

maximize
view([22 23]), light, material dull
saveas(gcf, [figures_dir, 'induced_gamma_left.tif'])
view([-33 24]), light, material dull
saveas(gcf, [figures_dir, 'induced_gamma_right.tif'])
view([-172 0]), light, material dull
saveas(gcf, [figures_dir, 'induced_gamma_left_medial.tif'])
view([-145 0]), light, material dull
saveas(gcf, [figures_dir, 'induced_gamma_right_medial.tif'])

% colorbar
dum = load([projectdir, 'results/TFR_group.mat'],'L');
cfgp2=[];
cfgp2.colormap = cfgp.funcolormap;
cfgp2.zlim = cfgp.funcolorlim;
figure; ft_singleplotTFR(cfgp2, dum.L);
saveas(gcf, [figures_dir, 'induced_gamma_colorbar.svg'])

%% decoding accuracy
addpath([projectdir, 'scripts/CanlabCore/CanlabCore/Visualization_functions/'])

x = load([projectdir, 'results/stat_decoding']); 
sprintf('DECODING ON MEG DATA')
sprintf('decoding accuracy attended, left %d percent (SD %d), p=%d; right %d percent (SD %d), p=%d', round(100*x.accuracy(1,1).mean), round(100*x.accuracy(1,1).std,1), x.stat(1,1).prob, round(100*x.accuracy(1,2).mean), round(100*x.accuracy(1,2).std,1), x.stat(1,2).prob)
sprintf('decoding accuracy unattended, left %d percent (SD %d), p=%d; right %d percent (SD %d), p=%d', round(100*x.accuracy(2,2).mean), round(100*x.accuracy(2,2).std,1), x.stat(2,2).prob, round(100*x.accuracy(2,1).mean), round(100*x.accuracy(2,1).std,1), x.stat(2,1).prob)

y = load([projectdir, 'results/stat_decoding_eye']); 
sprintf('DECODING ON EYE DATA')
sprintf('decoding accuracy attended, left %d percent (SD %d), p=%d; right %d percent (SD %d), p=%d', round(100*y.accuracy(1,1).mean), round(100*y.accuracy(1,1).std,1), y.stat(1,1).prob, round(100*y.accuracy(1,2).mean), round(100*y.accuracy(1,2).std,1), y.stat(1,2).prob)
sprintf('decoding accuracy unattended, left %d percent (SD %d), p=%d; right %d percent (SD %d), p=%d', round(100*y.accuracy(2,2).mean), round(100*y.accuracy(2,2).std,1), y.stat(2,2).prob, round(100*y.accuracy(2,1).mean), round(100*y.accuracy(2,1).std,1), y.stat(2,1).prob)

cmap = brewermap(2,'RdBu');
figure; violinplot(100*x.accuracy(1,1).all, 'facecolor', cmap(1,:), 'medc', [], 'x', 1, 'pointsize', 1)
hold on, violinplot(100*y.accuracy(1,1).all, 'facecolor', cmap(2,:), 'medc', [], 'x', 1, 'pointsize', 1)
violinplot(100*x.accuracy(1,2).all, 'facecolor', cmap(1,:), 'medc', [], 'x', 2, 'pointsize', 1)
violinplot(100*y.accuracy(1,2).all, 'facecolor', cmap(2,:), 'medc', [], 'x', 2, 'pointsize', 1)
ylim([45 80])
xlim([0 3])
saveas(gcf, [figures_dir, 'decoding_accuracy.eps'],'epsc')

%% Primals decoding

hemis=[1 2];
for subj=valid_subjects
  label = load([datadir, sprintf('sub%02d/meg%02d/sub%02d-meg%02d_cleandata.mat', subj, 2, subj, 2)]);
  label = label.data.label;
  for h=hemis
    filename = [results_dir, 'collapsephasebin/attended/', sprintf('sub%02d/sub%02d_decoding_hemi%d.mat', subj, subj, h)];
    tmp = load(filename, 'primal_P');
    primal_P{subj,h}.avg = tmp.primal_P;
    primal_P{subj,h}.label = label;
    primal_P{subj,h}.time = 0;
    primal_P{subj,h}.dimord = 'chan_time';
  end
end
cfgp=[];
cfgp.layout = 'CTF275_helmet.mat';
cfgp.colormap = flipud(brewermap(64, 'RdBu'));
cfgp.zlim = 'maxabs';
cfgp.numcontour = 0;
cfgp.gridscale = 250;
cfgp.marker = 'off';
cfgp.interpolation      = 'cubic'
sc = {'a', 'b'};
for h=hemis
  for subj=valid_subjects
    figure; maximize;
    ft_topoplotER(cfgp, primal_P{subj,h});
    saveas(gcf, [figures_dir, sprintf('primal_topo_new/primal_sub%02d%s.eps',subj, sc{h})], 'epsc')
  end
end



%% phasic modulation based on virtual channel
addpath('RainCloudPlots/tutorial_matlab/')
x=load([projectdir, 'results/stat_phasicmodulation_decoding']);
freqs=4:30;
side = {'left', 'right'};

for k=1:2
  [~, ix] = max(x.stat(k).stat);
  sprintf('attend %s max amp %s percent (SD = %s), p=%s at %d Hz', side{k}, num2str(round(mean(x.amp{k}(ix,:))*100,2)), num2str(round(std(x.amp{k}(ix,:))*100,2)), num2str(x.stat(k).uncorrected_p(ix)), freqs(ix))
end

cmap = (brewermap(2,'RdBu'));
figure
for k=1:2
  subplot(1,2,k)
  freqs=transpose(4:30);
  y = mean(mean(x.amprand{k},3),2);
  dy = mean(std(x.amprand{k},[],3),2);
  fill([freqs;flipud(freqs)],[y-dy;flipud(y+dy)],[0.8 0.8 0.8],'linestyle','none');
  line(freqs,y,'color', [0.5 0.5 0.5])
  hold on
  y = mean(x.amp{k},2);
  dy = std(x.amp{k},[],2);
  fill([freqs;flipud(freqs)],[y-dy;flipud(y+dy)],cmap(1,:),'linestyle','none');
  line(freqs,y,'color', 'r')
  xlim([4 30]), ylim([0 0.015])
  title(sprintf('attend %s', side{k}))
  saveas(gcf, [figures_dir, sprintf('phasic_modulation_virtualchan_%s.eps',side{k})],'eosc')
end

%% Phasic modulation of behavior
x = load([results_dir, 'stat_phasicmodulation_behavior']);
dum = load([projectdir, 'results/TFR_group.mat'],'L');

%%%%%%%%%%%%%
% statistic %
%%%%%%%%%%%%%
% attend left
x.stat{1}.mask = x.stat{1}.posclusterslabelmat==1; % visual inspection 
% shows that this cluster is mostly present in parietal regions, 9-15 Hz
f1 = find(x.stat{1}.time==9);
f2 = find(x.stat{1}.time==17);
params={'stat', 'std'};
for k=1:numel(params)
  y = x.stat{1}.(params{k}).*x.stat{1}.mask;
  y(y==0)=nan;
  s(k) = nanmean(nanmean(y(:,f1:f2)))*100;
end
sprintf('attend left: cluster 1: parietal regions, 9-17 Hz, clusterstat = %d, mean = %d , SD = %d percent, p = %d',100*x.stat{1}.posclusters(1).clusterstat, s(1), s(2), x.stat{1}.posclusters(1).prob)

% attend right
x.stat{2}.mask = x.stat{2}.posclusterslabelmat==1; % visual inspection 
% shows that this cluster is mostly present in parietal regions, 9-15 Hz
f1 = find(x.stat{1}.time==9);
f2 = find(x.stat{1}.time==17);
params={'stat', 'std'};
for k=1:numel(params)
  y = x.stat{2}.(params{k}).*x.stat{2}.mask;
  y(y==0)=nan;
  s(k) = nanmean(nanmean(y(:,f1:f2)))*100;
end
sprintf('attend right: cluster 1: parietal regions, 9-17 Hz, clusterstat = %d, mean = %d , SD = %d percent, p = %d',100*x.stat{2}.posclusters(1).clusterstat, s(1), s(2), x.stat{2}.posclusters(1).prob)


%%%%%%%%%%%%%%%%%%
% color settings %
%%%%%%%%%%%%%%%%%%
c = brewermap(10, 'Set3');
n=64;
cmap = (brewermap(64, 'RdBu'));
cmap = flipud(cmap(1:n/2,:));
cmap2 = (brewermap(2,'RdBu'));

%%%%%%%%%%%%%%%%%%%%%%%
% soueceplot settings %
%%%%%%%%%%%%%%%%%%%%%%%
cfgp=[];
cfgp.funparameter = 'stat';
cfgp.funcolormap = cmap;
cfgp.method = 'surface';
cfgp.colorbar = 'no';
cfgp.funcolorlim = [0.004 0.012];
cfgp.maskstyle = 'colormix';
cfgp.facecolor = 'skin';
cfgp.maskparameter = cfgp.funparameter;

%%%%%%%%%%%%%%%
% attend left %
%%%%%%%%%%%%%%%
k=1;
s = x.stat{k};
t=8;
tx = s.time==t;
s.stat = s.stat(:,tx);
s.time = t;
ft_sourceplot(cfgp, s)
view([0 90])
maximize
material dull
saveas(gcf, [figures_dir, 'modulation_behavior_left_8hz.tif'])

ix = match_str(s.label,  {'L_8_B05_06'});
figure; 
y=[2e-3 23e-3];
hold on, plot(4:30, x.stat{k}.randmean(ix,:), 'k')
plot(4:30, x.stat{k}.stat(ix,:), 'color', cmap2(1,:))
% first plot 8 Hz
hold on, for ii=1:10
plot(30+[1 5], [squeeze(x.ampr{k}(ix,tx,ii)), squeeze(x.amp{k}(ix,tx,ii))], '--o', 'color', c(ii,:))
end
plot(30+[1 5], [mean(x.ampr{k}(ix,tx,:)), mean(x.amp{k}(ix,tx,:))], '-ok', 'linewidth', 2)
% also plot 17 Hz
tx2 = x.stat{k}.time==17;
for ii=1:10
plot(30+[6 10], [squeeze(x.ampr{k}(ix,tx2,ii)), squeeze(x.amp{k}(ix,tx2,ii))], '--o', 'color', c(ii,:))
end
plot(30+[6 10], [mean(x.ampr{k}(ix,tx2,:)), mean(x.amp{k}(ix, tx2,:))], '-ok', 'linewidth',2 )
xlim([4 41]), ylim(y)
maximize
saveas(gcf, [figures_dir, 'modulation_behavior_left_L8_06.eps'],'epsc')

s = x.stat{k};
t=11;
tx = s.time==t;
s.stat = s.stat(:,tx);
s.time = t;
ft_sourceplot(cfgp, s)
view([0 90])
maximize
material dull
saveas(gcf, [figures_dir, 'modulation_behavior_left_11hz.tif'])

ix = match_str(s.label,  {'L_2_B05_06'});
figure; 
y=[2e-3 23e-3];
hold on, plot(4:30, x.stat{k}.randmean(ix,:), 'k')
plot(4:30, x.stat{k}.stat(ix,:), 'color', cmap2(1,:))
hold on, for ii=1:10
plot(30+[1 5], [squeeze(x.ampr{k}(ix,tx,ii)), squeeze(x.amp{k}(ix,tx,ii))], '--o', 'color', c(ii,:))
end
plot(30+[1 5], [mean(x.ampr{k}(ix,tx,:)), mean(x.amp{k}(ix,tx,:))], '-ok', 'linewidth', 2)
xlim([4 36]), ylim(y)
maximize
saveas(gcf, [figures_dir, 'modulation_behavior_left_L2_06.eps'],'epsc')

s = x.stat{k};
t=15;
tx = s.time==t;
s.stat = s.stat(:,tx);
s.time = t;
ft_sourceplot(cfgp, s)
view([0 60])
maximize
material dull
saveas(gcf, [figures_dir, 'modulation_behavior_left_15hz.tif'])

ix = match_str(s.label,  {'R_7_B05_01'});
figure; 
y=[0e-3 24e-3];
hold on, plot(4:30, x.stat{k}.randmean(ix,:), 'k')
plot(4:30, x.stat{k}.stat(ix,:), 'color', cmap2(1,:))
hold on, for ii=1:10
plot(30+[1 5], [squeeze(x.ampr{k}(ix,tx,ii)), squeeze(x.amp{k}(ix,tx,ii))], '--o', 'color', c(ii,:))
end
plot(30+[1 5], [mean(x.ampr{k}(ix,tx,:)), mean(x.amp{k}(ix,tx,:))], '-ok', 'linewidth', 2)
xlim([4 36]), ylim(y)
maximize
saveas(gcf, [figures_dir, 'modulation_behavior_left_R7_01.eps'],'epsc')

%%%%%%%%%%%%%%%%
% attend right %
%%%%%%%%%%%%%%%%
k=2;
s = x.stat{k};
t=13;
tx = s.time==t;
s.stat = s.stat(:,tx);
s.time = t;
ft_sourceplot(cfgp, s)
maximize
material dull
saveas(gcf, [figures_dir, 'modulation_behavior_right_13hz.tif'])

ix = match_str(s.label,  {'L_8_B05_02'});
figure; 
y=[3e-3 22e-3];
hold on, plot(4:30, x.stat{k}.randmean(ix,:), 'k')
plot(4:30, x.stat{k}.stat(ix,:), 'color', cmap2(1,:))
hold on, for ii=1:10
plot(30+[1 5], [squeeze(x.ampr{k}(ix,tx,ii)), squeeze(x.amp{k}(ix,tx,ii))], '--o', 'color', c(ii,:))
end
plot(30+[1 5], [mean(x.ampr{k}(ix,tx,:)), mean(x.amp{k}(ix,tx,:))], '-ok', 'linewidth', 2)
xlim([4 36]), ylim(y)
maximize
saveas(gcf, [figures_dir, 'modulation_behavior_right_L8_02.eps'],'epsc')

ix = match_str(s.label,  {'L_5_B05_03'});
figure; 
y=[0e-3 19e-3];
hold on, plot(4:30, x.stat{k}.randmean(ix,:), 'k')
plot(4:30, x.stat{k}.stat(ix,:), 'color', cmap2(1,:))
hold on, for ii=1:10
plot(30+[1 5], [squeeze(x.ampr{k}(ix,tx,ii)), squeeze(x.amp{k}(ix,tx,ii))], '--o', 'color', c(ii,:))
end
plot(30+[1 5], [mean(x.ampr{k}(ix,tx,:)), mean(x.amp{k}(ix,tx,:))], '-ok', 'linewidth', 2)
xlim([4 36]), ylim(y)
maximize
saveas(gcf, [figures_dir, 'modulation_behavior_right_L5_03.eps'],'epsc')

cfgp2=[];
cfgp2.colormap = cfgp.funcolormap;
cfgp2.zlim = cfgp.funcolorlim;
figure; ft_singleplotTFR(cfgp2, dum.L);
saveas(gcf, [figures_dir, 'modulation_behavior_colorbar.eps'],'epsc')

%% Modulation decoding - parcel level
x=load([projectdir, 'results/stat_phasicmodulation_decoding_parc']);
addpath([projectdir, 'scripts/CanlabCore/CanlabCore/Visualization_functions/'])

%%%%%%%%%%%%%
% statistic %
%%%%%%%%%%%%%
% attend left
x.stat{1}.mask = x.stat{1}.posclusterslabelmat==1; % visual inspection 
% shows that this concerns right FEF
f1 = find(x.stat{1}.time==17);
f2 = find(x.stat{1}.time==12);
params={'stat', 'std'};
for k=1:numel(params)
  y = x.stat{1}.(params{k}).*x.stat{1}.mask;
  y(y==0)=nan;
  s(k) = nanmean(nanmean(y(:,f1:f2)))*100;
end
sprintf('attend left: cluster 1: right FEF, 20 Hz, clusterstat = %d, mean = %d , SD = %d percent, p = %d', x.stat{1}.posclusters(1).clusterstat, s(1), s(2), x.stat{1}.posclusters(1).prob)

% attend right
x.stat{2}.mask = x.stat{2}.posclusterslabelmat==1 %| x.stat{2}.posclusterslabelmat==2; % visual inspection 
% shows that this cluster is mostly present in parietal regions, 9-15 Hz
f1 = find(x.stat{1}.time==8);
f2 = find(x.stat{1}.time==12);
params={'stat', 'std'};
for k=1:numel(params)
  y = x.stat{2}.(params{k}).*x.stat{2}.mask;
  y(y==0)=nan;
  s(k) = nanmean(nanmean(y(:,f1:f2)))*100;
end
sprintf('attend right: cluster 1: left parietal cortex and FEF, 8-12 Hz, clusterstat = %d, mean = %d , SD = %d percent, p = %d', 100*x.stat{2}.posclusters(1).clusterstat, s(1), s(2), x.stat{2}.posclusters(1).prob)


% color settings
c=brewermap(10,'Set3');
dum = load([projectdir, 'results/TFR_group.mat'],'L');
n=64;
cmap = (brewermap(n, 'RdBu'));
cmap = flipud(cmap(1:n/2,:));
cmap2 = (brewermap(2,'RdBu'));

% sourceplot settings
cfgp=[];
cfgp.funparameter = 'stat';
cfgp.funcolormap = cmap;
cfgp.method = 'surface';
cfgp.colorbar = 'no';
cfgp.maskstyle = 'colormix';
cfgp.facecolor = 'skin';
cfgp.maskparameter = cfgp.funparameter;
cfgp.funcolorlim = [0 0.014];

%%%%%%%%%%%%%%%
% attend left %
%%%%%%%%%%%%%%%
k=1;
s = x.stat{k};
t = 4;
tx = s.time==t;
s.time = t;
s.stat = s.stat(:,tx);

% sourceplot
ft_sourceplot(cfgp, s)
view([0 90])
maximize
material dull
saveas(gcf, [figures_dir, 'modulation_parc_left_4hz.tif'])

% line plot
figure;
ix = match_str(x.stat{k}.label, 'R_8_B05_06');
hold on, plot(4:20, x.stat{k}.randmean(ix,:), 'k')
plot(4:20, x.stat{k}.stat(ix,:), 'color', cmap2(1,:))
for ii=1:10
plot(20+[1 4], [squeeze(x.ampr{k}(ix,tx,ii)), squeeze(x.amp{k}(ix,tx,ii))], '--o', 'color', c(ii,:))
end
plot(20+[1 4], [mean(x.ampr{k}(ix,tx,:)), mean(x.amp{k}(ix,tx,:))], '-ok', 'linewidth', 2)
xlim([4 25]); ylim([5e-3 28e-3])
maximize
saveas(gcf, [figures_dir, 'modulation_parc_left_R8_06.eps'],'epsc')

s = x.stat{k};
t = 9;
tx = s.time==t;
s.time = t;
s.stat = s.stat(:,tx);

% sourceplot
ft_sourceplot(cfgp, s)
view([0 90])
maximize
material dull
saveas(gcf, [figures_dir, 'modulation_parc_left_9hz.tif'])

% line plot
figure;
ix = match_str(x.stat{k}.label, 'R_19_B05_12');
hold on, plot(4:20, x.stat{k}.randmean(ix,:), 'k')
plot(4:20, x.stat{k}.stat(ix,:), 'color', cmap2(1,:))
for ii=1:10
plot(20+[1 4], [squeeze(x.ampr{k}(ix,tx,ii)), squeeze(x.amp{k}(ix,tx,ii))], '--o', 'color', c(ii,:))
end
plot(20+[1 4], [mean(x.ampr{k}(ix,tx,:)), mean(x.amp{k}(ix,tx,:))], '-ok', 'linewidth', 2)
xlim([4 25]); ylim([2e-3 20e-3])
maximize
saveas(gcf, [figures_dir, 'modulation_parc_left_R19_12.eps'],'epsc')

s = x.stat{k};
t = 20;
tx = s.time==t;
s.time = t;
s.stat = s.stat(:,tx);

% sourceplot
ft_sourceplot(cfgp, s)
view([0 90])
maximize
material dull
saveas(gcf, [figures_dir, 'modulation_parc_left_20hz.tif'])

% line plot
figure;
ix = match_str(x.stat{k}.label, 'R_8_B05_02');
hold on, plot(4:20, x.stat{k}.randmean(ix,:), 'k')
plot(4:20, x.stat{k}.stat(ix,:), 'color', cmap2(1,:))
for ii=1:10
plot(20+[1 4], [squeeze(x.ampr{k}(ix,tx,ii)), squeeze(x.amp{k}(ix,tx,ii))], '--o', 'color', c(ii,:))
end
plot(20+[1 4], [mean(x.ampr{k}(ix,tx,:)), mean(x.amp{k}(ix,tx,:))], '-ok', 'linewidth', 2)
xlim([4 25]); ylim([2e-3 20e-3])
maximize
saveas(gcf, [figures_dir, 'modulation_parc_left_R8_02.eps'],'epsc')

% rose plot
usefreq = {4:6, 8:11, 16:20};
parc = {'R_8_B05_06', 'R_19_B05_12', 'R_8_B05_02'};
for jj=1:3
  figure;
  clear theta
  for subj=1:10
    ii=0;
    for f=usefreq{jj}
      ii=ii+1;
      [~,~,theta(subj, ii)] = analysis_phase_rtdiff(subj,k,f,parc{jj});
    end
    subplot(2,5,subj); hold on; rose(theta(subj,:)); 
  end
end

%%%%%%%%%%%%%%%%
% attend right %
%%%%%%%%%%%%%%%%
% #1
k=2;
s = x.stat{k};
t = 4;
tx = s.time==t;
s.time = t;
s.stat = s.stat(:,tx);

% sourceplot
ft_sourceplot(cfgp, s)
view([0 90])
maximize
material dull
saveas(gcf, [figures_dir, 'modulation_parc_right_4hz.tif'])

% line plot
figure;
ix = match_str(x.stat{k}.label, 'L_7_B05_05');
hold on, plot(4:20, x.stat{k}.randmean(ix,:), 'k')
plot(4:20, x.stat{k}.stat(ix,:), 'color', cmap2(1,:))
for ii=1:10
plot(20+[1 4], [squeeze(x.ampr{k}(ix,tx,ii)), squeeze(x.amp{k}(ix,tx,ii))], '--o', 'color', c(ii,:))
end
plot(20+[1 4], [mean(x.ampr{k}(ix,tx,:)), mean(x.amp{k}(ix,tx,:))], '-ok', 'linewidth', 2)
xlim([4 25]); ylim([4e-3 20e-3])
maximize
saveas(gcf, [figures_dir, 'modulation_parc_right_L7_05.eps'],'epsc')

% #2
s = x.stat{k};
t = 10;
tx = s.time==t;
s.time = t;
s.stat = s.stat(:,tx);

% sourceplot
ft_sourceplot(cfgp, s)
view([0 90])
maximize
material dull
saveas(gcf, [figures_dir, 'modulation_parc_right_10hz.tif'])

% line plot #1
figure;
ix = match_str(x.stat{k}.label, 'L_5_B05_02');
hold on, plot(4:20, x.stat{k}.randmean(ix,:), 'k')
plot(4:20, x.stat{k}.stat(ix,:), 'color', cmap2(1,:))
for ii=1:10
plot(20+[1 4], [squeeze(x.ampr{k}(ix,tx,ii)), squeeze(x.amp{k}(ix,tx,ii))], '--o', 'color', c(ii,:))
end
plot(20+[1 4], [mean(x.ampr{k}(ix,tx,:)), mean(x.amp{k}(ix,tx,:))], '-ok', 'linewidth', 2)
xlim([4 25]); ylim([2e-3 18e-3])
maximize
saveas(gcf, [figures_dir, 'modulation_parc_right_L5_02.eps'],'epsc')


% % line plot #2
% figure;
% ix = match_str(x.stat{k}.label, 'L_8_B05_06');
% hold on, plot(4:20, x.stat{k}.randmean(ix,:), 'k')
% plot(4:20, x.stat{k}.stat(ix,:), 'color', cmap2(1,:))
% for ii=1:10
% plot(20+[1 4], [squeeze(x.ampr{k}(ix,tx,ii)), squeeze(x.amp{k}(ix,tx,ii))], '--o', 'color', c(ii,:))
% end
% plot(20+[1 4], [mean(x.ampr{k}(ix,tx,:)), mean(x.amp{k}(ix,tx,:))], '-ok', 'linewidth', 2)
% xlim([4 25]); ylim([3e-3 16e-3])
% maximize
% saveas(gcf, [figures_dir, 'modulation_parc_right_L8_06.eps'],'epsc')

s = x.stat{k};
t = 13;
tx = s.time==t;
s.time = t;
s.stat = s.stat(:,tx);

% sourceplot
ft_sourceplot(cfgp, s)
view([0 90])
maximize
material dull
saveas(gcf, [figures_dir, 'modulation_parc_right_13hz.tif'])

% line plot #1
figure;
ix = match_str(x.stat{k}.label, 'L_8_B05_02');
hold on, plot(4:20, x.stat{k}.randmean(ix,:), 'k')
plot(4:20, x.stat{k}.stat(ix,:), 'color', cmap2(1,:))
for ii=1:10
plot(20+[1 4], [squeeze(x.ampr{k}(ix,tx,ii)), squeeze(x.amp{k}(ix,tx,ii))], '--o', 'color', c(ii,:))
end
plot(20+[1 4], [mean(x.ampr{k}(ix,tx,:)), mean(x.amp{k}(ix,tx,:))], '-ok', 'linewidth', 2)
xlim([4 25]); ylim([4e-3 20e-3])
maximize
saveas(gcf, [figures_dir, 'modulation_parc_right_L8_02.eps'],'epsc')

% colorbar
cfgp2=[];
cfgp2.colormap = cfgp.funcolormap;
cfgp2.zlim = cfgp.funcolorlim;
figure; ft_singleplotTFR(cfgp2, dum.L);
saveas(gcf, [figures_dir, 'modulation_parc_colorbar.eps'],'epsc')

%% Suplemental : phase
datainfo;
load atlas_subparc374_8k.mat
exclude_label = match_str(atlas.parcellationlabel, {'L_???_01', 'L_MEDIAL.WALL_01', 'R_???_01', 'R_MEDIAL.WALL_01'});
selparc = setdiff(1:numel(atlas.parcellationlabel),exclude_label); % hard coded exclusion of midline and ???

for k=1:numel(useparc)
  whichparc{k} = find(contains(atlas.parcellationlabel(selparc), useparc{k}));
end
whichparc = unique(cat(1,whichparc{:}));

c=brewermap(10,'Set3');
%%%%%%%%%%%%%%%
% attend left %
%%%%%%%%%%%%%%%
parc = {'R_8_B05_06','R_8_B05_02','R_8_B05_03'};
% parc = atlas.parcellationlabel(selparc(whichparc(19:21)));
clear theta
for jj=1:numel(parc)
for subj=1:10
ii=0;
for f=4:6
ii=ii+1;
[~,~,theta1(subj, ii,jj)] = analysis_phase_rtdiff(subj,1,f,parc{jj});
end
end
end

parc = {'R_19_B05_12', 'R_7_B05_08', 'R_7_B05_09', 'R_19_B05_09'};
% parc = atlas.parcellationlabel(selparc(whichparc(22:end)));
clear theta
for jj=1:numel(parc)
for subj=1:10
ii=0;
for f=8:10
ii=ii+1;
[~,~,theta2(subj, ii,jj)] = analysis_phase_rtdiff(subj,1,f,parc{jj});
end
end
end

parc = {'R_8_B05_02', 'R_8_B05_03', 'R_8_B05_06'};
% parc = atlas.parcellationlabel(selparc(whichparc(19:21)));
for jj=1:numel(parc)
for subj=1:10
ii=0;
for f=16:20
ii=ii+1;
[~,~,theta3(subj, ii,jj)] = analysis_phase_rtdiff(subj,1,f,parc{jj});
end
end
end

%%%%%%%%%%%%%%%%
% attend right %
%%%%%%%%%%%%%%%%
parc = {'L_7_B05_05', 'L_7_B05_01', 'L_7_B05_02', 'L_7_B05_09', 'L_7_B05_08', 'L_19_B05_09'};
% parc = atlas.parcellationlabel(selparc(whichparc(4:18)));
clear theta
for jj=1:numel(parc)
for subj=1:10
ii=0;
for f=4:6
ii=ii+1;
[~,~,theta4(subj, ii,jj)] = analysis_phase_rtdiff(subj,2,f,parc{jj});
end
end
end

parc = {'L_5_B05_02','L_2_B05_08'};
% parc = atlas.parcellationlabel(selparc(whichparc(4:18)));
clear theta
for jj=1:numel(parc)
for subj=1:10
ii=0;
for f=8:12
ii=ii+1;
[~,~,theta5(subj, ii,jj)] = analysis_phase_rtdiff(subj,2,f,parc{jj});
end
end
end

parc = {'L_8_B05_02', 'L_8_B05_03', 'L_8_B05_06'};
% parc = atlas.parcellationlabel(selparc(whichparc(1:3)));
for jj=1:numel(parc)
for subj=1:10
ii=0;
for f=10:18
ii=ii+1;
[~,~,theta6(subj, ii,jj)] = analysis_phase_rtdiff(subj,2,f,parc{jj});
end
end
end

N=18;
figure; for ii=1:10
subplot(2,5,ii); polarhistogram(reshape(theta1(ii,:,:),[],1),N, 'FaceColor', c(ii,:));
end
maximize
figure; for ii=1:10
subplot(2,5,ii); polarhistogram(reshape(theta2(ii,:,:),[],1),N, 'FaceColor', c(ii,:));
end
maximize
figure; for ii=1:10
subplot(2,5,ii); polarhistogram(reshape(theta3(ii,:,:),[],1),N, 'FaceColor', c(ii,:));
end
maximize
figure; for ii=1:10
subplot(2,5,ii); polarhistogram(reshape(theta4(ii,:,:),[],1),N, 'FaceColor', c(ii,:));
end
maximize
figure; for ii=1:10
subplot(2,5,ii); polarhistogram(reshape(theta5(ii,:,:),[],1),N, 'FaceColor', c(ii,:));
end
maximize
figure; for ii=1:10
subplot(2,5,ii); polarhistogram(reshape(theta6(ii,:,:),[],1),N, 'FaceColor', c(ii,:));
end
maximize
