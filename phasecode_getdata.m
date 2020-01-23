function [dat, dat_append] = phasecode_loaddata(subj, varargin)

doparc = ft_getopt(varargin, 'doparc', false);

datainfo;

if doparc
  load([projectdir, sprintf('results/tlck/sub%02d_sourceparc.mat', subj)])
end

cfg=[];
cnt=1;
for ses=subjects(subj).validsessions
    filename = [datadir, sprintf('sub%02d/meg%02d/sub%02d-meg%02d_cleandata.mat', subj, ses, subj, ses)];
    dat{cnt} = load(filename, 'data');
    dat{cnt} = removefields(dat{cnt}.data, 'elec');
    if doparc
        for k=1:numel(dat{cnt}.trial)
            dat{cnt}.trial{k} = source_parc{1}.F * dat{cnt}.trial{k};
        end
        dat{cnt}.label = source_parc{cnt}.label;
    end
    cnt=cnt+1;
end
cfg.appenddim = 'rpt';
dat_append = ft_appenddata(cfg, dat{:});
