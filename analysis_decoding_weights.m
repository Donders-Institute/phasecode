% this script is intended to visualize the weights from the decoding
% analysis.

datainfo;
% choose variables
subj=4;
hemi = 1;
f = 10;
contrast = 'attended';

%% load data
filename = [projectdir 'results/collapsephasebin/'];
switch contrast
  case 'congruent'
    filename = [filename, 'congruent/'];
  case 'attended'
    filename = [filename, 'attended/'];
  case 'unattended'
    filename = [filename, 'unattended/'];
end
filename = [filename, sprintf('sub%02d/sub%02d_decoding', subj, subj)];

if strcmp(contrast, 'attended') || strcmp(contrast, 'unattended')
  filename = [filename, sprintf('_hemi%d', hemi)];
end
filename = [filename, sprintf('_f%d', f)];
load(filename)


cfg=[];
cnt=1;
for ses=subjects(subj).validsessions
  filename = [datadir, sprintf('sub%02d/meg%02d/sub%02d-meg%02d_cleandata.mat', subj, ses, subj, ses)];
  tmp{cnt} = load(filename, 'data');
  tmp{cnt} = removefields(tmp{cnt}.data, 'elec');
  cnt=cnt+1;
end
data = ft_appenddata([], tmp{:});

%% create FT structure
dat.time = 1:18;
dat.label = data.label;
dat.dimord = 'chan_time';
dat.avg = squeeze(mean(mean(primal_P,4),2));
% observed that the topographies over folds and over permutations are
% qualitatively very similar. Over bins there is more variability. No
% interpretable topo though.

% plot
cfgp=[];
cfgp.layout = 'CTF275_helmet.mat';
cfgp.parameter = 'avg';
figure; ft_topoplotER(cfgp, dat)

%% Project to source level (2D)
x = load([projectdir, sprintf('results/tlck/sub%02d_sourceparc.mat', subj)]);
F = x.source_parc{hemi}.F;

load atlas_subparc374_8k.mat
load cortex_inflated_shifted.mat
atlas.pos=ctx.pos;
exclude_label = match_str(atlas.parcellationlabel, {'L_???_01', 'L_MEDIAL.WALL_01', 'R_???_01', 'R_MEDIAL.WALL_01'});
selparc = setdiff(1:numel(atlas.parcellationlabel),exclude_label); % hard coded exclusion of midline and ???

source=dat;
source.brainordinate = atlas;
source.label = atlas.parcellationlabel;
source.avg = zeros(374,18);
source.avg(selparc,:) = F*dat.avg;

cfgp2=rmfield(cfgp,'parameter');
cfgp2.funparameter = 'avg';
ft_sourcemovie(cfgp2, source);

%% cosinefit on svm weights
cfg
cfg=[];
cfg.cosinefit.statistic = 'complex';
s = statfun_cosinefit(cfg, source.avg, centerphase);

source2=source;
source2.time=0;
source2.avg = abs(s.stat);

cfgp2=cfgp;
cfgp2.method = 'surface';
cfgp2.funparameter = 'avg';
ft_sourceplot(cfgp2, source2);

%% Contrast with those of random permutations

filename = [projectdir 'results/collapsephasebin/'];
switch contrast
case 'congruent'
filename = [filename, 'congruent/'];
case 'attended'
filename = [filename, 'attended/'];
case 'unattended'
filename = [filename, 'unattended/'];
end
filename = [filename, sprintf('sub%02d/primalp/sub%02d_decoding_rand', subj, subj)];
if strcmp(contrast, 'attended') || strcmp(contrast, 'unattended')
filename = [filename, sprintf('_hemi%d', hemi)];
end
filename = [filename, sprintf('_f%d', f)];

for k=1:10
x{k} = load([filename, sprintf('_%d',k )], 'primal_P');
x{k}=x{k}.primal_P;
end

x=cat(2,x{:});
x=(cat(5,x{:}));
x=squeeze(mean(mean(x,4),2));

x=reshape(x,270,[]);
x=F*x;
x=reshape(x, [370,18,100]);

for k=1:100
tmp = statfun_cosinefit(cfg, x(:,:,k), centerphase);
srand(:,k) = abs(tmp.stat);
end

source3 = source2;
source3.avg(selparc,:) = (source3.avg(selparc,:)-mean(srand,2))./std(srand,[],2);
ft_sourceplot(cfgp2, source3);



%% Project to source level (3D)
x=load([projectdir, sprintf('results/freq/sub%02d_dics_gamma', subj)])
Fdics=x.dicsfilter{hemi};
Fdics = cat(1,Fdics{:});

load standard_sourcemodel3d6mm.mat
sourcemodel.avg = zeros(size(sourcemodel.inside));
sourcemodel.avg(sourcemodel.inside,1) = mean(Fdics*dat.avg,2);
