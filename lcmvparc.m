function [source_parc] = lcmvparc(data, tlck, F, L, atlas)

tmp     = rmfield(data, {'elec' 'grad','hdr','trialinfo', 'sampleinfo'});
selparc = 1:numel(atlas.parcellationlabel);
exclude_label = match_str(atlas.parcellationlabel, {'L_???_01', 'L_MEDIAL.WALL_01', 'R_???_01', 'R_MEDIAL.WALL_01'});
selparc = setdiff(1:numel(atlas.parcellationlabel),exclude_label); % hard coded exclusion of midline and ???

source_parc.label = atlas.parcellationlabel(selparc);
source_parc.time  = tlck.time;
source_parc.F     = cell(numel(source_parc.label),1);
source_parc.L     = cell(numel(source_parc.label),1);
source_parc.avg   = zeros(numel(selparc),numel(source_parc.time));
source_parc.dimord = 'chan_time';

% prepare the cfg for pca
cfg                       = [];
cfg.method                = 'pca';

for k = 1:numel(selparc)
    tmpF = F(atlas.parcellation==selparc(k),:);
    tmp.avg = tmpF*squeeze(mean(data.trial,1));
    for l=1:size(tmpF,1)
        tmplabel{l} = sprintf('loc%03d', l);
    end
    tmp.label = tmplabel;
    clear tmplabel
    cfg.comment = sprintf('Select the spatial filters for all nodes belonging to parcel %s. Create the source time courses for the nodes in this parcel by multiplying the spatial filters with the channel level data. Give this as input to ft_componentanalysis.', atlas.parcellationlabel{selparc(k)});
    tmpcomp   = ft_componentanalysis(cfg, tmp);

    tmpL = L(:,atlas.parcellation==selparc(k));
    
    source_parc.F{k}     = tmpcomp.unmixing*tmpF;
    source_parc.L{k}     = tmpL*tmpcomp.unmixing';
    source_parc.avg(k,:) = source_parc.F{k}(1,:)*squeeze(mean(tlck.trial,1));
    
    source_parc.F{k} = source_parc.F{k}(1,:);
end


source_parc.F = cat(1, source_parc.F{:});
