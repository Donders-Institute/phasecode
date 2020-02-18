function startup

if ispc
else
  addpath('/project/3011085.06/scripts/');
  addpath('/project/3011085.06/scripts/cellfunction/');
  addpath('/project/3011085.06/data/');
  addpath('/project/3011085.06/results/');
  addpath('/project/3011085.06/scripts/CircStat');
  addpath('/project/3011085.06/scripts/circshift_columns/')
  addpath('/project/3011085.06/scripts/fieldtrip');
  addpath('/project/3011085.06/scripts/fieldtrip/qsub');
  addpath('/project/3011085.06/scripts/fieldtrip/external/dmlt/');
  addpath('/project/3011085.06/scripts/fieldtrip/external/dmlt/external/glmnet/')
  addpath('/project/3011085.06/scripts/fieldtrip/external/dmlt/')
  addpath('/project/3011085.06/scripts/fieldtrip/external/dmlt/external/')
  addpath('/project/3011085.06/scripts/fieldtrip/external/dmlt/external/svm/')
  addpath('/project/3011085.06/scripts/fieldtrip/external/brewermap');
  cd('/project/3011085.06/scripts/')
end
ft_defaults

end