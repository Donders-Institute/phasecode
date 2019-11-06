function [ampl,theta] = fit_sine_nobins(y,phases)
%from Zoefel et al (2019)
%% BZ: This code was provided by Nicholas Bland, Queensland University, 16/04/2019

%% y is trial outcome
X = ones(size(y,1),3);
X(:,2) = sin(phases);
X(:,3) = cos(phases);
theta = X\y;

ampl = sqrt(theta(2)^2 + theta(3)^2);