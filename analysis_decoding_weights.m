% this script is intended to visualize the weights from the decoding
% analysis.

datainfo;
% choose variables
subj=5;
hemis = [1 2];
freqs = 4:1:20;
contrast = 'attended';

%% load data
cfg=[];
cnt=1;
for ses=subjects(subj).validsessions
  filename = [datadir, sprintf('sub%02d/meg%02d/sub%02d-meg%02d_cleandata.mat', subj, ses, subj, ses)];
  tmp{cnt} = load(filename, 'data');
  tmp{cnt} = removefields(tmp{cnt}.data, 'elec');
  cnt=cnt+1;
end
data = ft_appenddata([], tmp{:});

filename = [projectdir 'results/collapsephasebin/'];
switch contrast
  case 'congruent'
    filename = [filename, 'congruent/'];
  case 'attended'
    filename = [filename, 'attended/'];
  case 'unattended'
    filename = [filename, 'unattended/'];
end
filename_orig = [filename, sprintf('sub%02d/primalp/sub%02d_decoding', subj, subj)];
ii=0;
for f = freqs
  ii=ii+1;
  for h=hemis
    
if strcmp(contrast, 'attended') || strcmp(contrast, 'unattended')
  filename = [filename_orig, sprintf('_hemi%d', h)];
end
filename = [filename, sprintf('_f%d', f)];
load(filename)


% create FT structure
% observed that the topographies over folds and over permutations are
% qualitatively very similar. Over bins there is more variability. No
% interpretable topo though.

dat{h}.time = freqs;
dat{h}.label = data.label;
dat{h}.dimord = 'chan_time';
dat{h}.avg(:,ii) = mean(mean(mean(primal_P,4),3),2);

primP{h}(:,:,ii) = squeeze(mean(mean(primal_P,4),2));
  end
end

cfg=[];
cfg.cosinefit.statistic = 'complex';
nbin = 18;
centerphase = [0:2/nbin:(nbin-1)/(nbin/2)]*pi-pi;
for h=hemis
  [s1, s2, s3] = size(primP{h});
  s{h} = statfun_cosinefit(cfg, reshape(permute(primP{h}, [1 3 2]), [], s2), centerphase);
  dat{h}.stat = reshape(abs(s{h}.stat), [s1, s3]);
end


% plot
cfgp=[];
cfgp.layout = 'CTF275_helmet.mat';
cfgp.parameter = 'avg';
figure; ft_topoplotER(cfgp, dat{:})
cfgp.parameter = 'stat';
figure; ft_topoplotER(cfgp, dat{:})

%% Project to source level (2D)
x = load([projectdir, sprintf('results/tlck/sub%02d_sourceparc.mat', subj)]);

load atlas_subparc374_8k.mat
load cortex_inflated_shifted.mat
atlas.pos=ctx.pos;
exclude_label = match_str(atlas.parcellationlabel, {'L_???_01', 'L_MEDIAL.WALL_01', 'R_???_01', 'R_MEDIAL.WALL_01'});
selparc = setdiff(1:numel(atlas.parcellationlabel),exclude_label); % hard coded exclusion of midline and ???
for h=hemis
  F{h} = x.source_parc{h}.F;
  source{h}.brainordinate = atlas;
  source{h}.label = atlas.parcellationlabel;
  source{h}.avg = zeros(374,numel(freqs));
  source{h}.avg(selparc,:) = F{h}*dat{h}.avg;
  source{h}.stat(selparc,:) = F{h}*dat{h}.stat;
  source{h}.time = freqs;
end

cfgp2=rmfield(cfgp,'parameter');
cfgp2.funparameter = 'avg';
ft_sourcemovie(cfgp2, source);

