function startup

if ispc
    defaultdir = 'P:/';
else
    defaultdir = '/project/';
end
addpath([defaultdir, '3011085.06/scripts/']);
addpath([defaultdir, '3011085.06/scripts/cellfunction/']);
addpath([defaultdir, '3011085.06/data/']);
addpath([defaultdir, '3011085.06/results/']);
addpath([defaultdir, '3011085.06/scripts/CircStat']);
addpath([defaultdir, '3011085.06/scripts/circshift_columns/'])
addpath([defaultdir, '3011085.06/scripts/fieldtrip']);
addpath([defaultdir, '3011085.06/scripts/fieldtrip/qsub']);
addpath([defaultdir, '3011085.06/scripts/fieldtrip/external/dmlt/']);
addpath([defaultdir, '3011085.06/scripts/fieldtrip/external/dmlt/external/glmnet/'])
addpath([defaultdir, '3011085.06/scripts/fieldtrip/external/dmlt/'])
addpath([defaultdir, '3011085.06/scripts/fieldtrip/external/dmlt/external/'])
addpath([defaultdir, '3011085.06/scripts/fieldtrip/external/dmlt/external/svm/'])
addpath([defaultdir, '3011085.06/scripts/fieldtrip/external/brewermap']);
cd([defaultdir, '3011085.06/scripts/'])

ft_defaults

end