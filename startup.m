function startup

if ispc
else
  addpath('/project/3011085.02/phasecode/scripts/');
  addpath('/project/3011085.02/phasecode/scripts/cellfunction/');
  addpath('/project/3011085.02/phasecode/data/');
  addpath('/project/3011085.02/phasecode/results/');
  addpath('/project/3011085.02/phasecode/scripts/CircStat');
  addpath('/project/3011085.02/phasecode/scripts/circshift_columns/')
  addpath('/project/3011085.02/phasecode/scripts/fieldtrip');
  addpath('/project/3011085.02/phasecode/scripts/fieldtrip/qsub');
  addpath('/project/3011085.02/phasecode/scripts/fieldtrip/external/dmlt/');
  addpath('/project/3011085.02/phasecode/scripts/fieldtrip/external/dmlt/external/glmnet/')
  addpath('/project/3011085.02/phasecode/scripts/fieldtrip/external/dmlt/')
  addpath('/project/3011085.02/phasecode/scripts/fieldtrip/external/dmlt/external/')
  addpath('/project/3011085.02/phasecode/scripts/fieldtrip/external/dmlt/external/svm/')
  addpath('/project/3011085.02/phasecode/scripts/fieldtrip/external/brewermap');
  cd('/project/3011085.02/phasecode/scripts/')
end
ft_defaults

end