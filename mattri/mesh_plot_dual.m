function mesh_plot_dual(ax, mesh, bmarky)
% mesh_plot_dual(ax, mesh, bmarky)
%
% Plots a simplex-like mesh in red with the boundary edges in
% different color if option bmarky is 1

% This m-file is part of the mattri meshing package.
% Copyright (c) 2014, Mauro A Werder
% All rights reserved.
% Licensed under a BSD 2-Clause License, see LICENCE file

if nargin<3
    bmarky = 0;
end

nodes_vor = mesh.vor.nodes;
connect_edge_vor = mesh.vor.connect_edge;

axes(ax);
hold on;
vx = [nodes_vor(connect_edge_vor(:,1),1),nodes_vor(connect_edge_vor(:,2),1)].';
vy = [nodes_vor(connect_edge_vor(:,1),2),nodes_vor(connect_edge_vor(:,2),2)].';
h = plot(vx, vy, '-r'); %,nodes_vor(:,1),nodes_vor(:,2),'.r');
%set(h(1:end-1),'xliminclude','off','yliminclude','off')

if bmarky==1
    % edges with bmark=1
    bm1 = mesh.vor.bmark_edge==1;
    eg1a = mesh.vor.nodes(mesh.vor.connect_edge(bm1,1),:);
    eg1b = mesh.vor.nodes(mesh.vor.connect_edge(bm1,2),:);   
    eg1 = zeros(size(eg1a,1)*3,2)*NaN; % we need NaN inserted
                                       % between edges to stop
                                       % wrong edges being drawn
    eg1(1:3:end,:) = eg1a;
    eg1(2:3:end,:) = eg1b ;   
    plot(eg1(:,1),eg1(:,2), 'k')
    
    bm2 = mesh.vor.bmark_edge==2;
    eg1a = mesh.vor.nodes(mesh.vor.connect_edge(bm2,1),:);
    eg1b = mesh.vor.nodes(mesh.vor.connect_edge(bm2,2),:);   
    eg1 = zeros(size(eg1a,1)*3,2)*NaN;
    eg1(1:3:end,:) = eg1a;
    eg1(2:3:end,:) = eg1b ;   
    plot(eg1(:,1),eg1(:,2), 'g')
end


hold off;
axis equal;