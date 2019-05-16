function startup

if ispc
else
    addpath('/project/3011085.02/phasecode/scripts/');
    addpath('/project/3011085.02/phasecode/data/');
    addpath('/project/3011085.02/phasecode/results/');
    addpath('/home/common/matlab/fieldtrip');
    addpath('/home/common/matlab/fieldtrip/qsub');
    addpath('/home/common/matlab/fieldtrip/external/dmlt/');
    addpath('/home/common/matlab/fieldtrip/external/dmlt/external/glmnet/');
    cd('/project/3011085.02/phasecode/scripts/')
end
ft_defaults

end