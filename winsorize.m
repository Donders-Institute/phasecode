function [x] = winsorize(x, percent, dim)

% WINSORIZE winsorizes sample x in dimension dim setting the 
% percent percent observations in each of the tails
% to the value of the observation at threshold.
% Implemented according to Wilcox Robust estimation
% and hypothesis testing

if nargin<3,
  dim = find(size(x)>1, 1, 'first');
end

if nargin<2,
  percent = 0.2;
end

[thrlo, thrhi] = percthreshold(x, percent, dim);
n              = size(x, dim);
%g              = floor(percent*n);
repdimvec      = ones(1, ndims(x));
repdimvec(dim) = n;

%repdimvec(dim) = g;
%x(bsxfun(@lt, x, thrlo)) = repmat(thrlo, repdimvec);
%x(bsxfun(@gt, x, thrhi)) = repmat(thrhi, repdimvec);
x = max(x, repmat(thrlo, repdimvec));
x = min(x, repmat(thrhi, repdimvec));

