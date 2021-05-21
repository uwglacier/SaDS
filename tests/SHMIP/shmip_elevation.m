function z = shmip_elevation(xy)
% SHMIP_ELEVATION calculates the surface elevation using SHMIP geometry
%
% z = shmip_elevation(xy)
%
% Calculates elevation as
% z = 6*(sqrt(xy(:, 1) + 5e3) - sqrt(5e3));
%
% See also shmip_melt

z = 6*(sqrt(xy(:, 1) + 5e3) - sqrt(5e3));
end