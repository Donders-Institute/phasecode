function [s, transform] = shift_geometry(s, xrot, yrot, zrot, xshiftm, yshiftm, zshiftm)
% shifts the hemispheres of a source object such that they are next to each
% other. pos field is required in input (s). x,y,z are only to be specified
% if a different shift is required.
% Copyright (C) 202 Mats van Es, Donders Institute. m.vanes@donders.ru.nl

if ~exist('xrot', 'var'), xrot = 0; end
if ~exist('yrot', 'var'), yrot = 0; end
if ~exist('zrot', 'var'), zrot = pi/2; end


l = size(s.pos,1);
xrange = round(max(s.pos(:,1)) + abs(min(s.pos(:,1))));

for k=1:2
if ~exist('xshiftm', 'var'), xshift = (-1)^k * xrange; else xshift=xshiftm; end
if ~exist('yshiftm', 'var'), yshift = 0; else yshift=yshiftm; end
if ~exist('zshiftm', 'var'), zshift = 0; else zshift=zshiftm; end

ztmp = zrot*(-1)^k;    

cx = cos(xrot);
sx = sin(xrot);
cy = cos(yrot);
sy = sin(yrot);
cz = cos(ztmp);
sz = sin(ztmp);

% from http://air.bmap.ucla.edu/AIR5/homogenous.html
transform{k} = [(cz*cy+sz*sx*sy),  (sz*cy-cz*sx*sy), (cx*sy), xshift
                 (-sz*cx),           (cz*cx),        sx,      yshift
             (cz*sx*cy-cz*sy), (-cz*sx*cy-sz*sy), (cx*cy),    zshift
                     0                 0              0         0];
end
                 % s = ft_transform_geometry(transform, s);
if isfield(s, 'hemisphere')
  g1 = s.hemisphere==1;
  g2 = s.hemisphere==2;
  s.pos(g1,:) = ft_warp_apply(transform{1}, s.pos(g1,:));
  s.pos(g2,:) = ft_warp_apply(transform{2}, s.pos(g2,:));
else
  s.pos(1:l/2,:) = ft_warp_apply(transform{1}, s.pos(1:l/2,:));
  s.pos(l/2+1:end,:) = ft_warp_apply(transform{2}, s.pos(l/2+1:end,:));
end
