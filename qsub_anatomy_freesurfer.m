function qsub_anatomy_freesurfer(subj, script_number)


datainfo;
if script_number==1
    % freesurfer script 1
    shell_script = [projectdir, 'scripts/anatomy_freesurfer.sh'];
elseif script_number==2
    shell_script = [projectdir, 'scripts/anatomy_freesurfer2.sh'];
elseif script_number==3
    shell_script = [projectdir, 'scripts/anatomy_postfreesurferscript.sh'];
end

mri_dir = [projectdir, sprintf('results/anatomy/sub%02d/', subj)];
preproc_dir = 'preproc';

% create the string that is executed in the linux terminal
command = [shell_script, ' ', mri_dir, ' ', preproc_dir];

% call the script
system(command);


