% This m-file is part of the mattri meshing package.
% Copyright (c) 2014, Mauro A Werder
% All rights reserved.
% Licensed under a BSD 2-Clause License, see LICENCE file

addpath('..')
% make a bounding polygon
boundary_xy = [0,0;
               1,0;
               0.5,0.5;
               1,1;
               0,1;];  % do not close the loop!

% boundary marks>0 on edge:
bmark = [1;2;2;1;3];  % just a mark which is given to the nodes on the boundary
bmark_edge = [1;2;1;1;1];  % just a mark which is given to the edges on the boundary

maxarea = 0.01;       % maximal area of triangles

triareafn = @(xy) maxarea*(xy(:,1).^2 + xy(:,2).^2+0.01);
%triareafn = @(xy) maxarea*(xy(:,1).^2 + xy(:,2).^2+0.01);

[mesh, boundary_inds] = make_mesh(boundary_xy, bmark, bmark_edge, maxarea, triareafn);

figure;
mesh_plot_tri(gca, mesh.tri, 1, 1)
title('A mesh refined near the origin')
