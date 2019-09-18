 function [s, cfg] = statfun_cosinefit(cfg, dat, design)

% STATFUN_xxx is a function for computing a statistic for the relation
% between biological data and a design vector containing trial
% classifications or another independent variable
%
% This function is called by STATISTICS_MONTECARLO, where you can specify
% cfg.statistic = 'xxx' which will be evaluated as statfun_xxx.
%
% The external interface of this function has to be
%   [s] = statfun_xxx(cfg, dat, design);
% where
%   dat    contains the biological data, Nvoxels x Nreplications
%   design contains the independent variable,  1 x Nreplications
%
% Additional settings can be passed through to this function using
% the cfg structure.
%
% STATFUN_COSINEFIT fits a cosine to the data dat. the independent
% variable design should contain angular values, bounded by -pi and pi.
%
% the output s is a structure containing the statistic, as specified by
%  cfg.cosinefit.statistic. this can either be the amplitude (default) of the fit, angle   (giving the preferred angle), complex (giving both angle and amplitude in a complex number), or fit (giving the percentage of explained variance).
% additional fields in the output=structure are:
%  s.r      percentage of explained variance.
%  s.offset the DC-component of the fit.
%
% Additional cfg-options are:
%  cfg.cosinefit.phi = [] (default) or angular value between -pi and pi, estimates the amplitude of a cosine at a given angle.

% Copyright (C) 2006 Jan-Mathijs Schoffelen
%
% $Log: statfun_cosinefit.m,v $
% Revision 1.2  2009/10/12 14:56:08  jansch
% some changes to speed up
%
% Revision 1.1  2008/09/23 07:31:14  roboos
% Moved all statfuns and trialfuns to their own directories, where they will be easier to find for the end-user. Also updated fieldtripdefs accordingly.
%
% Revision 1.2  2007/08/01 12:23:18  jansch
% optimise functions for matrices as input
%
% Revision 1.1  2006/04/13 11:15:00  jansch
% first implementation for fieldtrip
%

if ~isfield(cfg, 'cosinefit'),           cfg.cosinefit           = [];          end
if ~isfield(cfg.cosinefit, 'phi'),       cfg.cosinefit.phi       = [];          end
if ~isfield(cfg.cosinefit, 'statistic'), cfg.cosinefit.statistic = 'amplitude'; end

%---create and check independent variable
if size(design, 2) ~= size(dat,2), error('design is incompatible with input-data'); end
if size(design, 1) ~= size(dat,1) && size(design,1)==1,
  quickflag = 1;
  repdim    = [size(dat,1) 1];
else
  quickflag = 0;
  repdim    = [1 1];
end  

%---estimate offset
offset = mean(dat,2);

%---do cosinefitting based on a least-square fit
y = dat - repmat(offset, [1 size(dat,2)]);

n  = size(design,2);
if quickflag,
  S  = y*sin(design(:));
  C  = y*cos(design(:));
  dS = repmat(sum(sin(2.*design)), repdim);
  dC = repmat(sum(cos(2.*design)), repdim);
else
  S  = sum(y.*sin(design),2);
  C  = sum(y.*cos(design),2);
  dS = sum(sin(2.*design),2);
  dC = sum(cos(2.*design),2);
end

if isempty(cfg.cosinefit.phi),
  nom   = S.*n + S.*dC - C.*dS; 
  denom = C.*n - S.*dS - C.*dC;
  b(denom > 0,1) = atan(nom(denom > 0)./denom(denom > 0));
  b(denom < 0,1) = atan(nom(denom < 0)./denom(denom < 0)) + pi;
  b(denom ==0,1) = 0.5 .* pi;
else
  b = ones(size(design,1),1).*cfg.cosinefit.phi;
end
  
A    = 2 .* (C.*cos(b)+S.*sin(b)) ./ (dS.*cos(2.*b)+dC.*sin(2.*b)+n);
sse  = sum( (y - repmat(A, [1 n]).*cos(repmat(design,repdim)-repmat(b, [1 n]))).^2, 2); %sum of squared errors
sTot = sum(y.^2, 2); 
sA   = sqrt( (sse./(n-3)) ./ sum( (design - repmat(mean(design,2), [1 n])).^2, 2) );  %standard deviation of estimator
r    = 1 - sse ./ sTot; %proportion of explained variance
phi  = b; 

%---create output-structure
if strcmp(cfg.cosinefit.statistic, 'complex'),
  s.stat = A.*exp(1i.*phi);
elseif strcmp(cfg.cosinefit.statistic, 'amplitude'),
  s.stat = A;
elseif strcmp(cfg.cosinefit.statistic, 'angle'),
  s.stat = phi;
elseif strcmp(cfg.cosinefit.statistic, 'fit'),
  s.stat = r;
elseif strcmp(cfg.cosinefit.statistic, 'tstat'),
  s.stat = A./sA;
%elseif strcmp(cfg.cosinefit.statistic, 'fstat'),
%  s.stat = sse./(sTot./(n-3));
elseif strcmp(cfg.cosinefit.statistic, 'ratio'),
  s.stat = A./offset;
end
s.r      = r;
%s.s      = sA;
s.offset = offset;
