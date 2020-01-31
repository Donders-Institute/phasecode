
datainfo;
addpath([projectdir, 'scripts/kldiv/']);
contrast = 'attended';
model = '3d';
if ~exist('method', 'var'), method = 'cosinefit'; end

for subj=1:10
    subj
    if strcmp(model, '2d')
        doparc = true;
    else
      doparc = false;
    end
[~, data_all] = phasecode_getdata(subj, 'doparc', doparc);
    
    freqs=4:1:30;
    cnt=0;
    for f = 1:numel(freqs)
        cnt=cnt+1;
        for hemi = 1:2
            [data, idx] = phasecode_select_condition(data_all, contrast, hemi);
            filename = [projectdir, 'results/phase/', sprintf('sub%02d_phase_%d', subj, freqs(f))];
            load(filename)
            phasebin = phasebin{mod(hemi,2)+1}; % for attention left, take right hemisphere phase.
            phase = phase{mod(hemi,2)+1};
            % end
            phase(~idx,:)=[];
            
            % only select validly cued and correct trials
            cfg=[];
            cfg.trials = (data.trialinfo(:,7)==1) & (data.trialinfo(:,1)==data.trialinfo(:,4));
            data = ft_selectdata(cfg, data);
            rmtrials = find(~cfg.trials);
            phase(rmtrials,:,:)=[];
            
            % take RT
            rt = data.trialinfo(:,6);
            
            % take phase at 1 sec
            x = nearest(time,1);
            phi = phase(:, x);
            
            % bin
            nbins=50;
            bins = [0:2/nbins:2*(nbins-1)/nbins]*pi;
            rtbin = zeros(nbins,1);
            nperm = 1000;
            rtbinshuff = zeros(nbins, nperm);
            for k=1:nbins
                centerphase = bins(k);
                lim = [centerphase-pi/2 centerphase+pi/2];
                lim(lim<0) = lim(lim<0) + 2*pi;
                lim(lim>2*pi) = lim(lim>2*pi) - 2*pi;
                if lim(1)>lim(2)
                    sel = phi>lim(1) | phi<lim(2);
                else
                    sel = phi>lim(1) & phi<lim(2);
                end
                rtbin(k,1) = mean(rt(sel));
                
                for ii=1:nperm
                    rtshuff = rt(randperm(numel(rt)));
                    rtbinshuff(k,ii) = mean(rtshuff(sel));
                end
            end
            
            switch method
                case 'cosinefit'
                    cfg=[];
                    cfg.cosinefit.statistic = 'complex';
                    zx = numel(rtbin);
                    centerphase = [0:2/zx:(zx-1)/(zx/2)]*pi-pi;
                    s = statfun_cosinefit(cfg, rtbin', centerphase);
                    amp(cnt, hemi) = abs(s.stat);
                    ang(cnt, hemi) = angle(s.stat);
                    for k=1:nperm
                        srand = statfun_cosinefit(cfg, rtbinshuff(:,k)', centerphase);
                        amprand(k, cnt, hemi) = abs(srand.stat);
                        angrand(k, cnt, hemi) = angle(srand.stat);
                    end
                case 'kldiv'
                    % compute kl-divergence (deviation from uniformity)
                    uniform = ones(nbins,1)/nbins;
                    rtbin_distr = rtbin./sum(rtbin);
                    kl(cnt,hemi) = kldiv([], uniform, rtbin_distr);
                    for k=1:nperm
                        rtbinshuff_distr = rtbinshuff(:,k)./sum(rtbinshuff(:,k));
                        klshuff(k) = kldiv([], uniform, rtbinshuff_distr);
                    end
                    kl_u(cnt,hemi) = mean(klshuff);
                    kl_std(cnt,hemi) = std(klshuff);
            end
        end
    end
    
    switch method
        case 'cosinefit'
            all(subj).amp = amp;
            all(subj).ang = ang;
            all(subj).amprand = amprand;
            all(subj).angrand = angrand;
            keep all subj contrast model subjects projectdir datadir method
        case 'kldiv'
            all(subj).kl=kl;
            all(subj).kl_u=kl_u;
            all(subj).kl_std=kl_std;
            keep all subj contrast model subjects projectdir datadir method
    end
end

%% channel/parcel level
datainfo;
addpath([projectdir, 'scripts/kldiv/']);
contrast = 'attended';
method = 'cosinefit';
subj=4;
model='2d';
if strcmp(model, '2d')
    doparc = true;
    load atlas_subparc374_8k.mat
else
  doparc = false;
end
[~, data] = phasecode_getdata(subj, 'doparc', doparc);

cfg=[];
cfg.trials = (data.trialinfo(:,7)==1) & (data.trialinfo(:,1)==data.trialinfo(:,4));
data = ft_selectdata(cfg, data);

cnt=0;
freqs=4:1:20;
for f=freqs
    cnt=cnt+1;
    fs = 200;
    % mtmconvol
    numcycles = 2;
    cfg=[];
    cfg.method = 'mtmfft';
    cfg.taper = 'hanning';
    cfg.foi = f;
    cfg.pad = 4;
    cfg.output = 'fourier';
    cfg.t_ftimwin = numcycles*1./cfg.foi;
    cfg2.latency = [1-(numcycles/2*1/f) 1+(numcycles/2*1/f)];
    cfg.keeptrials = 'yes';
    freq = ft_freqanalysis(cfg, ft_selectdata(cfg2, data));
    %     time = cfg.toi;
    %     tcnt(cnt) = nearest(freq.time,1);
    phase{cnt} = angle(squeeze(freq.fourierspctrm))+pi;
end

rt = data.trialinfo(:,6);
phi=cat(3,phase{:});
hemis= [1 2];
nbins=50;
nperm = 100;
bins = [0:2/nbins:2*(nbins-1)/nbins]*pi;
nfreq = numel(freqs);

for f=1:nfreq
    f
    rtbin{f} = zeros(nbins, numel(hemis), size(phi,2));
    rtbinshuff{f} = zeros(nbins, numel(hemis), size(phi,2), nperm);
    for hemi=hemis
        trlidx = find(data.trialinfo(:,1)==hemi);
        tmpphi = phi(trlidx,:,f);
        tmprt = rt(trlidx);

        for k=1:nbins
            centerphase = bins(k);
            lim = [centerphase-pi/2 centerphase+pi/2];
            lim(lim<0) = lim(lim<0) + 2*pi;
            lim(lim>2*pi) = lim(lim>2*pi) - 2*pi;
            if lim(1)>lim(2)
                sel = tmpphi>lim(1) | tmpphi<lim(2);
            else
                sel = tmpphi>lim(1) & tmpphi<lim(2);
            end
            if isempty(sel), break
            end
            for ch = 1:size(sel,2)
                rtbin{f}(k,hemi,ch) = mean(tmprt(sel(:,ch)));
            end
            
            for ii=1:nperm
                rtshuff = tmprt(randperm(numel(tmprt)));
                for ch = 1:size(sel,2)
                    rtbinshuff{f}(k,hemi,ch,ii) = mean(rtshuff(sel(:,ch)));
                end
            end
            
        end
    end
end

switch method
    case 'cosinefit'
        % compute cosine fit
        cfg=[];
        cfg.cosinefit.statistic = 'complex';
        zx = 50;
        centerphase = [0:2/zx:(zx-1)/(zx/2)]*pi-pi;
        
        amp = zeros(nfreq, numel(hemis), size(phi,2));
        ang = zeros(nfreq, numel(hemis), size(phi,2));
        amprand = zeros(nperm, nfreq, numel(hemis), size(phi,2));
        angrand = zeros(nperm, nfreq, numel(hemis), size(phi,2));
        
        for f=1:nfreq
            for h=hemis
                for ch=1:numel(data.label)
                    dat = rtbin{f}(:,h,ch);
                    s = statfun_cosinefit(cfg, dat', centerphase);
                    amp(f,h,ch) = abs(s.stat);
                    ang(f,h,ch) = angle(s.stat);
                    
                    tmp2 = squeeze(rtbinshuff{f}(:,h,ch,:));
                    for k=1:nperm
                        srand = statfun_cosinefit(cfg, tmp2(:,k)', centerphase);
                        amprand(k, f, h, ch) = abs(srand.stat);
                        angrand(k, f, h, ch) = angle(srand.stat);
                    end
                end
            end
        end
        
        y = (amp-squeeze(mean(amprand)))./squeeze(std(amprand));
        
    case 'kldiv'
        % compute kl-divergence (deviation from uniformity)
        uniform = ones(nbins,1)/nbins;
        for f=1:numel(rtbin)
            for h=hemis
                for ch=1:numel(data.label)
                    dat = rtbin{f}(:,h,ch);
                    tmp2 = squeeze(rtbinshuff{f}(:,h,ch,:));
                    rtbin_distr = dat/sum(dat);
                    kl(f,h,ch) = kldiv([], uniform, rtbin_distr);
                    for k=1:nperm
                        rtbinshuff_distr = tmp2(:,k)./sum(tmp2(:,k));
                        klshuff(k) = kldiv([], uniform, rtbinshuff_distr);
                    end
                    kl_u(f,h,ch) = mean(klshuff);
                    kl_std(f,h,ch) = std(klshuff);
                end
            end
        end
        
        y=(kl-kl_u)./kl_std;
end


stat{1} = rmfield(data, {'trial', 'time'});
stat{1}.dimord ='chan_time';
if strcmp(model, '2d')
    stat{1}.time=freqs;
    exclude_label = match_str(atlas.parcellationlabel, {'L_???_01', 'L_MEDIAL.WALL_01', 'R_???_01', 'R_MEDIAL.WALL_01'});
    selparc = setdiff(1:numel(atlas.parcellationlabel),exclude_label); % hard coded exclusion of midline and ???
    stat{1}.powspctrm = zeros(374,numel(freqs));
    stat{1}.powspctrm(selparc,:) = squeeze(y(:,1,:))';
    stat{1}.brainordinate = atlas;
    load cortex_inflated_shifted.mat
    stat{1}.brainordinate.pos = ctx.pos;
    stat{1}.label = atlas.parcellationlabel;
    stat{2}=stat{1}; stat{2}.powspctrm(selparc,:) = squeeze(y(:,2,:))';
else
    stat{1}.freq=freqs;
    stat{1}.powspctrm = squeeze(y(:,1,:))';
    stat{2}=stat{1}; stat{2}.powspctrm = squeeze(y(:,2,:))';
end


cfgp.layout = 'CTF275_helmet.mat';
cfgp.colormap = flipud(brewermap(64, 'RdBu'));
if strcmp(model,'2d')
    cfgp.funparameter = 'powspctrm';
    ft_sourcemovie(cfgp, stat{1})
else
    figure; ft_topoplotER(cfgp, stat{1}, stat{2})
end

save(sprintf('sub%02d_cosinefit_behavior_parc.mat', subj), 'stat', 'ang','amp', 'amprand', 'angrand','phase')

