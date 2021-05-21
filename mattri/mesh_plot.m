function mesh_plot(ax, mesh, labels)
% mesh_plot(ax, mesh, labels)
%
% Plots a mesh
% if labels == 1 then label the nodes

% This m-file is part of the mattri meshing package.
% Copyright (c) 2014, Mauro A Werder
% All rights reserved.
% Licensed under a BSD 2-Clause License, see LICENCE file

if nargin<3
    labels = 0;
end

mesh_plot_tri(ax, mesh.tri,1, labels);
hold on;
mesh_plot_simplex(ax, mesh.si,1);
hold off;
if mesh.scaled
    xlabel('x')
    ylabel('y')
else
    xlabel('x (m)')
    ylabel('y (m)')
end    