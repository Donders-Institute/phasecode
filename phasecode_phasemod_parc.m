subj=4;
datainfo;

tdir = [projectdir, sprintf('results/collapsephasebin/attended/sub%02d/parc/', subj)];

hemis = [1 2];
f = 10;
npermfiles = 2;
nparc = 370;

% load data
load atlas_subparc374_8k.mat
load cortex_inflated_shifted.mat
exclude_label = match_str(atlas.parcellationlabel, {'L_???_01', 'L_MEDIAL.WALL_01', 'R_???_01', 'R_MEDIAL.WALL_01'});
selparc = setdiff(1:numel(atlas.parcellationlabel),exclude_label); % hard coded exclusion of midline and ???


for h=hemis
  for k=1:nparc
    tmp1{k} = load([tdir, sprintf('sub%02d_decoding_hemi%d_f10_%d',subj,h, k)], 'accuracy');
    tmp1{k} = tmp1{k}.accuracy;
    
    for l=1:npermfiles
      tmp2{k,l} = load([tdir, sprintf('sub%02d_decoding_rand_hemi%d_f10_%d_%d',subj,h,k,l)], 'accuracy');
      tmp2{k,l} = tmp2{k,l}.accuracy;
    end
    tmp2{k,1} = cat(1,tmp2{k,:});
  end
  tmp2(:,2)=[];
  
  acc{h} = squeeze(mean(mean(cat(4, tmp1{:}),3),1));
  accr{h} = squeeze(mean(mean(cat(5,tmp2{:}),4),2));
  
  clear tmp1 tmp2
end


% cosinefit
nbin = size(acc{1},1);
centerphase = [0:2/nbin:(nbin-1)/(nbin/2)]*pi-pi;

cfg=[];
cfg.cosinefit.statistic = 'complex';
for h=hemis
  s{h} = statfun_cosinefit(cfg, acc{h}', centerphase);
  amp{h} = abs(s{h}(:).stat);
  ang{h} = angle(s{h}(:).stat);
  
  tmp = reshape(permute(accr{h}, [2,3,1]), nbin,[]);
  sr{h} = statfun_cosinefit(cfg, tmp', centerphase);
  ampr{h} = reshape(abs(sr{h}(:).stat), nparc, []);
end
  
% make fieldtrip structure

source=[];
source.brainordinate = atlas;
source.brainordinate.pos = ctx.pos;
source.label = atlas.parcellationlabel;
source.dimord = 'chan_time';
source.time = 1:numel(hemis);
source.stat = zeros(numel(source.label), numel(hemis));
for h=hemis
  source.stat(selparc,h) = (amp{h}-mean(ampr{h},2))./std(ampr{h},[],2);
end


% plot

cfgp = [];
cfgp.funparameter = 'stat';
ft_sourcemovie(cfgp, source);



  
  