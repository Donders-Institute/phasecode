

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
hold on, violinplot(100*x.accuracy(1,2).all, 'facecolor', cmap(1,:), 'medc', [], 'x', 2, 'pointsize', 1)
violinplot(100*y.accuracy(1,1).all, 'facecolor', cmap(2,:), 'medc', [], 'x', 1, 'pointsize', 1)
violinplot(100*y.accuracy(1,2).all, 'facecolor', cmap(2,:), 'medc', [], 'x', 2, 'pointsize', 1)
ylim([45 80])
xlim([0 3])
saveas(gcf, [figures_dir, 'decoding_accuracy.eps'],'epsc')

m1 = x.accuracyrand(1,1).mean;
s1 = x.accuracyrand(1,1).std;
m2 = x.accuracyrand(1,2).mean;
s2 = x.accuracyrand(1,2).std;
m3 = y.accuracyrand(1,1).mean;
s3 = y.accuracyrand(1,1).std;
m4 = y.accuracyrand(1,2).mean;
s4 = y.accuracyrand(1,2).std;
figure; shadedErrorBar(1:4, 100*[m1 m1 m2 m2], 100*[s1 s1 s2 s2]);ylim([45 80]); xlim([0 3])
saveas(gcf, [figures_dir, 'decoding_accuracy_rand.eps'],'epsc')
figure; shadedErrorBar(1:4, 100*[m3 m3 m4 m4], 100*[s3 s3 s4 s4]);ylim([45 80]); xlim([0 3])
saveas(gcf, [figures_dir, 'decoding_accuracy_rand_eye.eps'],'epsc')



%% phasic modulation based on virtual channel
addpath('RainCloudPlots/tutorial_matlab/')
x=load([projectdir, 'results/stat_phasicmodulation_decoding']);
freqs=4:30;
side = {'left', 'right'};

for k=1:2
  [~, ix(k)] = max(x.stat(k).stat);
  sprintf('attend %s max amp %s percent (SD = %s), p=%s at %d Hz', side{k}, num2str(round(mean(x.amp{k}(ix(k),:))*100,2)), num2str(round(std(x.amp{k}(ix(k),:))*100,2)), num2str(x.stat(k).uncorrected_p(ix(k))), freqs(ix(k)))
end

cmap = (brewermap(2,'RdBu'));
c = brewermap(10, 'dark2');
ix = [8 6];% 11 and 9 Hz
for k=1:2
  figure;
  subplot(1,2,1)
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
  
  subplot(1,2,2), hold on
  for ii=1:10
    plot([1 4], [squeeze(x.amprand{k}(ix(k),ii)), squeeze(x.amp{k}(ix(k),ii))], '--o', 'color', c(ii,:))
  end
  plot([1 4], [mean(x.amprand{k}(ix(k),:)), mean(x.amp{k}(ix(k),:))], '-ok', 'linewidth', 2)
  xlim([0 5]); ylim([0 0.02])
  maximize
%   saveas(gcf, [figures_dir, sprintf('phasic_modulation_virtualchan_%s.eps',side{k})],'epsc')
end

%% Phasic modulation of behavior
x = load([results_dir, 'stat_phasicmodulation_behavior']);
dum = load([projectdir, 'results/TFR_group.mat'],'L');
freqs = x.stat{1}.time;

% attend left
h=1;
x.stat{h}.mask = x.stat{h}.posclusterslabelmat==1; % visual inspection 
y = reshape(x.amp{h}, [],10);
y = mean(y(x.stat{h}.mask,:),1);
mu = mean(y);
sigma2 = std(y);
f1 = min(x.stat{h}.time(find(sum(x.stat{h}.posclusterslabelmat==1))));
f2 = max(x.stat{h}.time(find(sum(x.stat{h}.posclusterslabelmat==1))));
parc = x.stat{h}.label(find(sum(x.stat{h}.posclusterslabelmat==1,2)));
sprintf('attend left: cluster 1: %s, %d-%d Hz, clusterstat = %d, mean = %d , SD = %d percent, p = %d; BF = %d', parc{1}, f1, f2, x.stat{h}.posclusters(1).clusterstat, mu, sigma2, x.stat{h}.posclusters(1).prob, x.bf(h).bf10)

% attend right
h=2;
x.stat{h}.mask = x.stat{h}.posclusterslabelmat==1; % visual inspection 
y = reshape(x.amp{h}, [],10);
y = mean(y(x.stat{h}.mask,:),1);
mu = mean(y);
sigma2 = std(y);
f1 = min(x.stat{h}.time(find(sum(x.stat{h}.posclusterslabelmat==1))));
f2 = max(x.stat{h}.time(find(sum(x.stat{h}.posclusterslabelmat==1))));
parc = x.stat{h}.label(find(sum(x.stat{h}.posclusterslabelmat==1,2)));
sprintf('attend right: cluster 1: %s, %d-%d Hz, clusterstat = %d, mean = %d , SD = %d percent, p = %d; BF = %d', parc{1}, f1, f2, x.stat{h}.posclusters(1).clusterstat, mu, sigma2, x.stat{h}.posclusters(1).prob, x.bf(h).bf10)


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
Y = x.stat{k}.randmean';
dY = x.stat{k}.randstd';
Z = x.stat{k}.stat';
dZ = x.stat{k}.std';

s = x.stat{k};
t=8;
tx = s.time==t;
s.stat = s.mean(:,tx);
s.time = t;
ft_sourceplot(cfgp, s)
view([0 90])
maximize
material dull
saveas(gcf, [figures_dir, 'modulation_behavior_left_8hz.tif'])

ix = match_str(s.label,  {'L_8_B05_06'});
figure; hold on,
y=[2e-3 23e-3];
fill([freqs';flipud(freqs')],[Y(:,ix)-dY(:,ix);flipud(Y(:,ix)+dY(:,ix))],[0.8 0.8 0.8],'linestyle','none');
line(freqs,Y(:,ix),'color', [0.5 0.5 0.5])
fill([freqs';flipud(freqs')],[Z(:,ix)-dZ(:,ix);flipud(Z(:,ix)+dZ(:,ix))],cmap2(1,:),'linestyle','none');
line(freqs,Z(:,ix),'color', 'r')
% first plot 8 Hz
hold on, for ii=1:10
plot(freqs(end)+[1 5], [squeeze(x.ampr{k}(ix,tx,ii)), squeeze(x.amp{k}(ix,tx,ii))], '--o', 'color', c(ii,:))
end
plot(freqs(end)+[1 5], [mean(x.ampr{k}(ix,tx,:)), mean(x.amp{k}(ix,tx,:))], '-ok', 'linewidth', 2)
xlim([4 26]), ylim(y)
maximize
saveas(gcf, [figures_dir, 'modulation_behavior_left_L8_06.eps'],'epsc')

s = x.stat{k};
t=14;
tx = s.time==t;
s.stat = s.mean(:,tx);
s.time = t;
ft_sourceplot(cfgp, s)
view([0 90])
maximize
material dull
saveas(gcf, [figures_dir, 'modulation_behavior_left_14hz.tif'])

ix = match_str(s.label,  {'L_5_B05_01'});
figure; hold on
y=[2e-3 23e-3];
fill([freqs';flipud(freqs')],[Y(:,ix)-dY(:,ix);flipud(Y(:,ix)+dY(:,ix))],[0.8 0.8 0.8],'linestyle','none');
line(freqs,Y(:,ix),'color', [0.5 0.5 0.5])
fill([freqs';flipud(freqs')],[Z(:,ix)-dZ(:,ix);flipud(Z(:,ix)+dZ(:,ix))],cmap2(1,:),'linestyle','none');
line(freqs,Z(:,ix),'color', 'r')
hold on, for ii=1:10
plot(freqs(end)+[1 5], [squeeze(x.ampr{k}(ix,tx,ii)), squeeze(x.amp{k}(ix,tx,ii))], '--o', 'color', c(ii,:))
end
plot(freqs(end)+[1 5], [mean(x.ampr{k}(ix,tx,:)), mean(x.amp{k}(ix,tx,:))], '-ok', 'linewidth', 2)
xlim([4 26]), ylim(y)
maximize
saveas(gcf, [figures_dir, 'modulation_behavior_left_L5_01.eps'],'epsc')

s = x.stat{k};
t=12;
tx = s.time==t;
s.stat = s.mean(:,tx);
s.time = t;
ft_sourceplot(cfgp, s)
view([0 60])
maximize
material dull
saveas(gcf, [figures_dir, 'modulation_behavior_left_12hz.tif'])

ix = match_str(s.label,  {'R_7_B05_01'});
figure; hold on,
y=[2e-3 20e-3];
fill([freqs';flipud(freqs')],[Y(:,ix)-dY(:,ix);flipud(Y(:,ix)+dY(:,ix))],[0.8 0.8 0.8],'linestyle','none');
line(freqs,Y(:,ix),'color', [0.5 0.5 0.5])
fill([freqs';flipud(freqs')],[Z(:,ix)-dZ(:,ix);flipud(Z(:,ix)+dZ(:,ix))],cmap2(1,:),'linestyle','none');
line(freqs,Z(:,ix),'color', 'r')
hold on, for ii=1:10
plot(freqs(end)+[1 5], [squeeze(x.ampr{k}(ix,tx,ii)), squeeze(x.amp{k}(ix,tx,ii))], '--o', 'color', c(ii,:))
end
plot(freqs(end)+[1 5], [mean(x.ampr{k}(ix,tx,:)), mean(x.amp{k}(ix,tx,:))], '-ok', 'linewidth', 2)
xlim([4 26]), ylim(y)
maximize
saveas(gcf, [figures_dir, 'modulation_behavior_left_R7_01.eps'],'epsc')

%%%%%%%%%%%%%%%%
% attend right %
%%%%%%%%%%%%%%%%
k=2;
Y = x.stat{k}.randmean';
dY = x.stat{k}.randstd';
Z = x.stat{k}.stat';
dZ = x.stat{k}.std';

s = x.stat{k};
t=13;
tx = s.time==t;
s.stat = s.mean(:,tx);
s.time = t;
ft_sourceplot(cfgp, s)
maximize
material dull
saveas(gcf, [figures_dir, 'modulation_behavior_right_13hz.tif'])

ix = match_str(s.label,  {'L_8_B05_02'});
figure; hold on,
y=[3e-3 22e-3];
fill([freqs';flipud(freqs')],[Y(:,ix)-dY(:,ix);flipud(Y(:,ix)+dY(:,ix))],[0.8 0.8 0.8],'linestyle','none');
line(freqs,Y(:,ix),'color', [0.5 0.5 0.5])
fill([freqs';flipud(freqs')],[Z(:,ix)-dZ(:,ix);flipud(Z(:,ix)+dZ(:,ix))],cmap2(1,:),'linestyle','none');
line(freqs,Z(:,ix),'color', 'r')
for ii=1:10
plot(freqs(end)+[1 5], [squeeze(x.ampr{k}(ix,tx,ii)), squeeze(x.amp{k}(ix,tx,ii))], '--o', 'color', c(ii,:))
end
plot(freqs(end)+[1 5], [mean(x.ampr{k}(ix,tx,:)), mean(x.amp{k}(ix,tx,:))], '-ok', 'linewidth', 2)
xlim([4 26]), ylim(y)
maximize
saveas(gcf, [figures_dir, 'modulation_behavior_right_L8_02.eps'],'epsc')

cfgp2=[];
cfgp2.colormap = cfgp.funcolormap;
cfgp2.zlim = cfgp.funcolorlim;
figure; ft_singleplotTFR(cfgp2, dum.L);
saveas(gcf, [figures_dir, 'modulation_behavior_colorbar.eps'],'epsc')

%% Modulation decoding - parcel level
x=load([projectdir, 'results/stat_phasicmodulation_decoding_parc']);
x.topostat = analysis_phasic_modulation_topography;
addpath([projectdir, 'scripts/CanlabCore/CanlabCore/Visualization_functions/'])

%%%%%%%%%%%%%
% statistic %
%%%%%%%%%%%%%
% attend left
h=1;
x.stat{h}.mask = x.stat{h}.posclusterslabelmat==1; % visual inspection 
% shows that this concerns right FEF
y = reshape(x.amp{h}, [],10);
y = mean(y(x.stat{h}.mask,:),1);
mu = mean(y);
sigma2 = std(y);
f1 = min(x.stat{h}.time(find(sum(x.stat{h}.posclusterslabelmat==1))));
f2 = max(x.stat{h}.time(find(sum(x.stat{h}.posclusterslabelmat==1))));
parc = x.stat{h}.label(find(sum(x.stat{h}.posclusterslabelmat==1,2)));
sprintf('attend left: cluster 1: %s, %d-%d Hz, clusterstat = %d, mean = %d , SD = %d percent, p = %d; BF = %d', parc{1}, f1, f2, x.stat{h}.posclusters(1).clusterstat, mu, sigma2, x.stat{h}.posclusters(1).prob, x.bf(h).bf10)

% attend right
h=2;
x.stat{h}.mask = x.stat{h}.posclusterslabelmat==1; % visual inspection 
% shows that this concerns right FEF
y = reshape(x.amp{h}, [],10);
y = mean(y(x.stat{h}.mask,:),1);
mu = mean(y);
sigma2 = std(y);
f1 = min(x.stat{h}.time(find(sum(x.stat{h}.posclusterslabelmat==1))));
f2 = max(x.stat{h}.time(find(sum(x.stat{h}.posclusterslabelmat==1))));
parc = x.stat{h}.label(find(sum(x.stat{h}.posclusterslabelmat==1,2)));
sprintf('attend right: cluster 1: %s, %d-%d Hz, clusterstat = %d, mean = %d , SD = %d percent, p = %d; BF = %d', parc{1}, f1, f2, x.stat{h}.posclusters(1).clusterstat, mu, sigma2, x.stat{h}.posclusters(1).prob, x.bf(h).bf10)


% color settings
c=brewermap(10,'Set3');
dum = load([projectdir, 'results/TFR_group.mat'],'L');
n=10;
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

% create mask showing the selected ROIs
load atlas_subparc374_8k.mat
for k=1:numel(useparc)
whichparc2{k} = find(contains(atlas.parcellationlabel, useparc{k}));
end
whichparc2 = unique(cat(1,whichparc2{:}));
x=ismember(atlas.parcellation, whichparc2);
figure; ft_plot_mesh(dum.brainordinate, 'contour', x, 'contourcolor', 'r', 'facecolor', 'skin')
lighting gouraud
camlight
material dull
view([0 90])
maximize
saveas(gcf, [figures_dir, 'methods_roi_mask.tif'])
%%%%%%%%%%%%%%%
% attend left %
%%%%%%%%%%%%%%%
k=1;
Y = x.stat{k}.randmean';
dY = x.stat{k}.randstd';
Z = x.stat{k}.stat';
dZ = x.stat{k}.std';

s = x.stat{k};
t = 4;
tx = s.time==t;
s.time = t;
s.stat = s.stat(:,tx);

% sourceplot
ft_sourceplot(cfgp, s2)
view([0 90])
maximize
material dull
saveas(gcf, [figures_dir, 'modulation_parc_left_4hz.tif'])

% line plot
figure;hold on
ix = match_str(x.stat{k}.label, 'R_8_B05_06');
fill([freqs';flipud(freqs')],[Y(:,ix)-dY(:,ix);flipud(Y(:,ix)+dY(:,ix))],[0.8 0.8 0.8],'linestyle','none');
line(freqs,Y(:,ix),'color', [0.5 0.5 0.5])
fill([freqs';flipud(freqs')],[Z(:,ix)-dZ(:,ix);flipud(Z(:,ix)+dZ(:,ix))],cmap2(1,:),'linestyle','none');
line(freqs,Z(:,ix),'color', 'r')
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
figure; hold on
ix = match_str(x.stat{k}.label, 'R_19_B05_12');
fill([freqs';flipud(freqs')],[Y(:,ix)-dY(:,ix);flipud(Y(:,ix)+dY(:,ix))],[0.8 0.8 0.8],'linestyle','none');
line(freqs,Y(:,ix),'color', [0.5 0.5 0.5])
fill([freqs';flipud(freqs')],[Z(:,ix)-dZ(:,ix);flipud(Z(:,ix)+dZ(:,ix))],cmap2(1,:),'linestyle','none');
line(freqs,Z(:,ix),'color', 'r')
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
figure;hold on
ix = match_str(x.stat{k}.label, 'R_8_B05_02');
fill([freqs';flipud(freqs')],[Y(:,ix)-dY(:,ix);flipud(Y(:,ix)+dY(:,ix))],[0.8 0.8 0.8],'linestyle','none');
line(freqs,Y(:,ix),'color', [0.5 0.5 0.5])
fill([freqs';flipud(freqs')],[Z(:,ix)-dZ(:,ix);flipud(Z(:,ix)+dZ(:,ix))],cmap2(1,:),'linestyle','none');
line(freqs,Z(:,ix),'color', 'r')
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
Y = x.stat{k}.randmean';
dY = x.stat{k}.randstd';
Z = x.stat{k}.stat';
dZ = x.stat{k}.std';

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
figure;hold on
ix = match_str(x.stat{k}.label, 'L_7_B05_05');
fill([freqs';flipud(freqs')],[Y(:,ix)-dY(:,ix);flipud(Y(:,ix)+dY(:,ix))],[0.8 0.8 0.8],'linestyle','none');
line(freqs,Y(:,ix),'color', [0.5 0.5 0.5])
fill([freqs';flipud(freqs')],[Z(:,ix)-dZ(:,ix);flipud(Z(:,ix)+dZ(:,ix))],cmap2(1,:),'linestyle','none');
line(freqs,Z(:,ix),'color', 'r')
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
figure;hold on
ix = match_str(x.stat{k}.label, 'L_5_B05_02');
fill([freqs';flipud(freqs')],[Y(:,ix)-dY(:,ix);flipud(Y(:,ix)+dY(:,ix))],[0.8 0.8 0.8],'linestyle','none');
line(freqs,Y(:,ix),'color', [0.5 0.5 0.5])
fill([freqs';flipud(freqs')],[Z(:,ix)-dZ(:,ix);flipud(Z(:,ix)+dZ(:,ix))],cmap2(1,:),'linestyle','none');
line(freqs,Z(:,ix),'color', 'r')
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
figure; hold on
ix = match_str(x.stat{k}.label, 'L_8_B05_02');
fill([freqs';flipud(freqs')],[Y(:,ix)-dY(:,ix);flipud(Y(:,ix)+dY(:,ix))],[0.8 0.8 0.8],'linestyle','none');
line(freqs,Y(:,ix),'color', [0.5 0.5 0.5])
fill([freqs';flipud(freqs')],[Z(:,ix)-dZ(:,ix);flipud(Z(:,ix)+dZ(:,ix))],cmap2(1,:),'linestyle','none');
line(freqs,Z(:,ix),'color', 'r')
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


%%%%%%%%%%%%%%
% WHOLEBRAIN %
clim = [0.015, 0.011, 0.011; 0.015, 0.011, 0.011];
for h=[1 2]
s=x.topostat{h};
s.time=0;
for k=1:3
s.stat = s.mean(:,k);
cfgp.funcolorlim = [0.006 clim(h,k)];
cfgp.funcolormap = cmap;
ft_sourceplot(cfgp, s)
view([0 90])
maximize
material dull
end
end

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
