

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
saveas(gcf, [figures_dir, 'TFR_topo_low.eps'])

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

cfgp=rmfield(cfgp, 'channel');
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
saveas(gcf, [figures_dir, 'TFR_topo_high.eps'])

cfg=[];
cfg.avgoverfreq = 'yes';
cfg.avgovertime = 'yes';
cfg.avgoverchan = 'yes';
cfg.latency = [0.4 1];
cfg.frequency = [50 65];
cfg.channel = chans;
H_avg = ft_selectdata(cfg, H);

sprintf('mean gamma power increase %d percent, STD %d percent', H_avg.powspctrm*100, H_avg.std*100)

%% source location gamma virtual channel
load atlas_subparc374_8k.mat
increase = [];
for subj=1:10
tmp = load([projectdir, sprintf('results/freq/sub%02d_dics_gamma', subj)], 'maxidx', 'powRatio');
loc(subj,:) = tmp.maxidx;
pow(subj,:) = tmp.powRatio.mean-1;
increase = [increase; 100*pow(subj,tmp.maxidx)];
end
sprintf('gamma power increase, mean left %d (SD %d percent)| right %d percent (SD %d percent)', round(mean(increase(:,1))), round(std(increase(:,1))), round(mean(increase(:,2))), round(std(increase(:,2))))

load standard_sourcemodel3d6mm.mat
load cortex_inflated_shifted.mat
P = sourcemodel;
P = shift_geometry(P);
% P.pos=ctx.pos;
P.avg = 100*nanmean(pow)';

virtualchanpos=[];
virtualchanpos.chanpos = P.pos(loc,:);
virtualchanpos.elecpos = P.pos(loc,:);
for k=1:20
  virtualchanpos.label{k}   = sprintf('%d', k);
end

cfgp=[];
cfgp.method = 'surface';
cmap = flipud(brewermap(128, 'RdBu'));
cmap(1:64,:) = [];
cfgp.funcolormap = cmap;
cfgp.funparameter = 'avg';
cfgp.funcolorlim = [0 80];
cfgp.surffile = 'surface_white_both_shifted.mat';
ft_sourceplot(cfgp, P); hold on
ft_plot_sens(virtualchanpos,'elecshape','sphere', 'elecsize', 0.5, 'facecolor', 'black')
view([0 52])
h=light;lighting gouraud
saveas(gcf, [figures_dir, 'induced_gamma.svg'])

dum = load([projectdir, 'results/TFR_group.mat'],'L')
cfgp2=[];
cfgp2.colormap = cfgp.funcolormap;
cfgp2.parameter = cfgp.funparameter;
cfgp2.colorlim = cfgp.funcolorlim;
figure; ft_singleplotTFR(cfgp2, dum.L);
saveas(gcf, [figures_dir, 'induced_gamma_colorbar.eps'])

%% decoding accuracy
addpath([projectdir, 'scripts/CanlabCore/CanlabCore/Visualization_functions/'])

x = load([projectdir, 'results/stat_decoding']); 
sprintf('DECODING ON MEG DATA')
sprintf('decoding accuracy attended, left %d percent (SD %d), p=%d; right %d percent (SD %d), p=%d', round(100*x.accuracy(1,1).mean), round(100*x.accuracy(1,1).std), x.stat(1,1).prob, round(100*x.accuracy(1,2).mean), round(100*x.accuracy(1,2).std), x.stat(1,2).prob)
sprintf('decoding accuracy unattended, left %d percent (SD %d), p=%d; right %d percent (SD %d), p=%d', round(100*x.accuracy(2,1).mean), round(100*x.accuracy(2,1).std), x.stat(2,1).prob, round(100*x.accuracy(2,2).mean), round(100*x.accuracy(2,2).std), x.stat(2,2).prob)

y = load([projectdir, 'results/stat_decoding_eye']); 
sprintf('DECODING ON EYE DATA')
sprintf('decoding accuracy attended, left %d percent (SD %d), p=%d; right %d percent (SD %d), p=%d', round(100*y.accuracy(1,1).mean), round(100*y.accuracy(1,1).std), y.stat(1,1).prob, round(100*y.accuracy(1,2).mean), round(100*y.accuracy(1,2).std), y.stat(1,2).prob)
sprintf('decoding accuracy unattended, left %d percent (SD %d), p=%d; right %d percent (SD %d), p=%d', round(100*y.accuracy(2,1).mean), round(100*y.accuracy(2,1).std), y.stat(2,1).prob, round(100*y.accuracy(2,2).mean), round(100*y.accuracy(2,2).std), y.stat(2,2).prob)

cmap = brewermap(2,'RdBu');
figure; violinplot(100*x.accuracy(1,1).all, 'facecolor', cmap(1,:), 'medc', [], 'x', 1, 'pointsize', 1)
hold on, violinplot(100*y.accuracy(1,1).all, 'facecolor', cmap(2,:), 'medc', [], 'x', 1, 'pointsize', 1)
violinplot(100*x.accuracy(1,2).all, 'facecolor', cmap(1,:), 'medc', [], 'x', 2, 'pointsize', 1)
violinplot(100*y.accuracy(1,2).all, 'facecolor', cmap(2,:), 'medc', [], 'x', 2, 'pointsize', 1)
ylim([45 80])
xlim([0 3])
saveas(gcf, [figures_dir, 'decoding_accuracy.eps'])

%% phasic modulation based on virtual channel
addpath('RainCloudPlots/tutorial_matlab/')
x=load([projectdir, 'results/stat_phasicmodulation_decoding']);
freqs=4:30;
leg = {'left', 'right'};

for k=1:2
[~, ix] = max(x.stat(k).stat);
sprintf('attend %s max amp %s percent (SD = %s), p=%s at %d Hz', leg{k}, num2str(round(mean(x.amp{k}(ix,:))*100,2)), num2str(round(std(x.amp{k}(ix,:))*100,2)), num2str(x.stat(k).uncorrected_p(ix)), freqs(ix))
end

figure;
cmap = (brewermap(2,'RdBu'));
for k=1:2
  subplot(1,2,k)
  hold on, shadedErrorBar(4:30, mean(mean(x.amprand{k},3),2), mean(std(x.amprand{k},[],3),2), [], 1)
  shadedErrorBar(4:30, mean(x.amp{k},2), std(x.amp{k},[],2), {'r', 'markerfacecolor', cmap(1,:)}, 1)
  title(sprintf('attend %s', leg{k}))
  xlim([4 30])
end
saveas(gcf, [figures_dir, 'phasic_modulation_virtualchan.eps'])

%% Phasic modulation of behavior
x = load([results_dir, 'stat_phasicmodulation_behavior']);
dum = load([projectdir, 'results/TFR_group.mat'],'L');

% color settings
c = brewermap(10, 'Set3');
n=64;
cmap = (brewermap(64, 'RdBu'));
cmap = flipud(cmap(1:n/2,:));
cmap2 = (brewermap(2,'RdBu'));

% soueceplot settings
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
saveas(gcf, [figures_dir, 'modulation_behavior_left_L8_06.eps'])

s = x.stat{k};
t=11;
tx = s.time==t;
s.stat = s.stat(:,tx);
s.time = t;
ft_sourceplot(cfgp, s)
view([0 90])
maximize
saveas(gcf, [figures_dir, 'modulation_behavior_left_11hz.tif'])

ix = match_str(s.label,  {'L_2_B05_06'});
figure; 
y=[2e-3 20e-3];
hold on, plot(4:30, x.stat{k}.randmean(ix,:), 'k')
plot(4:30, x.stat{k}.stat(ix,:), 'color', cmap2(1,:))
hold on, for ii=1:10
plot(30+[1 5], [squeeze(x.ampr{k}(ix,tx,ii)), squeeze(x.amp{k}(ix,tx,ii))], '--o', 'color', c(ii,:))
end
plot(30+[1 5], [mean(x.ampr{k}(ix,tx,:)), mean(x.amp{k}(ix,tx,:))], '-ok', 'linewidth', 2)
xlim([4 36]), ylim(y)
maximize
saveas(gcf, [figures_dir, 'modulation_behavior_left_L2_06.eps'])

s = x.stat{k};
t=15;
tx = s.time==t;
s.stat = s.stat(:,tx);
s.time = t;
ft_sourceplot(cfgp, s)
view([0 60])
maximize
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
saveas(gcf, [figures_dir, 'modulation_behavior_left_R7_01.eps'])

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
saveas(gcf, [figures_dir, 'modulation_behavior_right_L8_02.eps'])

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
saveas(gcf, [figures_dir, 'modulation_behavior_right_L5_03.eps'])

cfgp2=[];
cfgp2.colormap = cfgp.funcolormap;
cfgp2.zlim = cfgp.funcolorlim;
figure; ft_singleplotTFR(cfgp2, dum.L);
saveas(gcf, [figures_dir, 'modulation_behavior_colorbar.eps'])

%% Modulation decoding - parcel level
x=load([projectdir, 'results/stat_phasicmodulation_decoding_parc']);
addpath([projectdir, 'scripts/CanlabCore/CanlabCore/Visualization_functions/'])

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
tx = s.time==4;
s.time = t;
s.stat = s.stat(:,tx);

% sourceplot
ft_sourceplot(cfgp, s)
view([0 90])
maximize
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
saveas(gcf, [figures_dir, 'modulation_parc_left_R8_06.eps'])

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
xlim([4 25]); ylim([3e-3 27e-3])
maximize
saveas(gcf, [figures_dir, 'modulation_parc_left_L7_05.eps'])

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
saveas(gcf, [figures_dir, 'modulation_parc_right_10hz.tif'])

% line plot #1
figure;
ix = match_str(x.stat{k}.label, 'L_8_B05_06');
hold on, plot(4:20, x.stat{k}.randmean(ix,:), 'k')
plot(4:20, x.stat{k}.stat(ix,:), 'color', cmap2(1,:))
for ii=1:10
plot(20+[1 4], [squeeze(x.ampr{k}(ix,tx,ii)), squeeze(x.amp{k}(ix,tx,ii))], '--o', 'color', c(ii,:))
end
plot(20+[1 4], [mean(x.ampr{k}(ix,tx,:)), mean(x.amp{k}(ix,tx,:))], '-ok', 'linewidth', 2)
xlim([4 25]); ylim([3e-3 16e-3])
maximize
saveas(gcf, [figures_dir, 'modulation_parc_left_L8_06.eps'])

% line plot #2
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
saveas(gcf, [figures_dir, 'modulation_parc_left_L5_02.eps'])

% colorbar
cfgp2=[];
cfgp2.colormap = cfgp.funcolormap;
cfgp2.zlim = cfgp.funcolorlim;
figure; ft_singleplotTFR(cfgp2, dum.L);
saveas(gcf, [figures_dir, 'modulation_parc_colorbar.eps'])