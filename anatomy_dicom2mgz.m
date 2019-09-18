function mri = anatomy_dicom2mgz(subj)
% anatomy_dicom2mgz takes the the subject info data structure (or subject string as 'sXX') 
%   
%   Picks up the dicom files, reslices the image and creates a .mgz file (spm coordsyst)

datainfo;
filename = subjects(subj).mri;
mridir = [projectdir, sprintf('results/anatomy/sub%02d/', subj)];
mri = ft_read_mri(subjects(subj).mri);


% save the images in the mgz format
cfg             = [];
cfg.filename    = fullfile(mridir, 'mri.mgz');
cfg.filetype    = 'mgz';
cfg.parameter   = 'anatomy';
ft_volumewrite(cfg, mri);

end

