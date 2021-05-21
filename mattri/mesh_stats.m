function mesh_stats(mesh)
%  mesh_stats(mesh)
%
% Prints statistics about the mesh.
% plots a historgram of edge direction.

% This m-file is part of the mattri meshing package.
% Copyright (c) 2014, Mauro A Werder
% All rights reserved.
% Licensed under a BSD 2-Clause License, see LICENCE file

disp(' ')
disp('Triangular mesh')
disp('---------------')
disp(['Vertices: ', num2str(size(mesh.tri.nodes,1))])
disp(['Elements: ', num2str(size(mesh.tri.connect,1))])
disp(['Edges: ', num2str(size(mesh.tri.connect_edge,1))])

% figure out directions ignoring the edges at the boundary
bmark_e = mesh.tri.bmark_edge==0;
edge_coords1 =  mesh.tri.nodes(mesh.tri.connect_edge(bmark_e,1),:);
edge_coords2 =  mesh.tri.nodes(mesh.tri.connect_edge(bmark_e,2),:);
dir = edge_coords2-edge_coords1;
dir_angle = atan(dir(:,2)./dir(:,1));

figure
hist(dir_angle,30);
title('Distribution of edge orientation angles')
xlabel('Angle (rad)')
ylabel('number')


disp(' ')
disp('Dual mesh')
disp('---------')
disp(['Vertices: ', num2str(size(mesh.vor.nodes,1))])
disp(['Cells: ', num2str(size(mesh.vor.connect,1))])
disp(['Edges: ', num2str(size(mesh.vor.connect_edge,1))])