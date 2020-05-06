function [rt_opt_mean, rt_sub_mean, optimal, rt_opt, rt_sub] = analysis_phase_rtdiff(subj, hemi, f, parc)
datainfo
freqs = 4:20;
fix = find(freqs==f);

load atlas_subparc374_8k.mat
exclude_label = match_str(atlas.parcellationlabel, {'L_???_01', 'L_MEDIAL.WALL_01', 'R_???_01', 'R_MEDIAL.WALL_01'});
selparc = setdiff(1:numel(atlas.parcellationlabel),exclude_label); % hard coded exclusion of midline and ???

% load single trial RT and phase
filename = [projectdir, 'results/modulation/', sprintf('sub%02d_cosinefit_behavior_parc.mat', subj)];
load(filename, 'condition', 'phi', 'rt')

ix1 = match_str(atlas.parcellationlabel(selparc), parc);
idx = condition==hemi;
phi = phi(idx, ix1, fix)-pi; % phi was manually added in a previous analysis. get back to [-pi pi] range.
rt = rt(idx);

% load optimal phase
x=load([projectdir, 'results/stat_phasicmodulation_decoding_parc'], 'ang');
ix2 = match_str(atlas.parcellationlabel, parc);
optimal = x.ang{hemi}(ix2,fix,subj);

% get the distance from optimal phase
dist = circ_dist(phi, optimal);
idx_opt = dist<=pi/2;
idx_sub = dist>pi/2;

% get the average RT for (sub)optimal phase
rt_opt = rt(idx_opt);
rt_sub = rt(idx_sub);

rt_opt_mean = mean(rt_opt);
rt_sub_mean = mean(rt_sub);



