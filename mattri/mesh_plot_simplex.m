function mesh_plot_simplex(ax,meshsi, bmarky)
% mesh_plot_simplex(ax,meshsi, bmarky)
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


nodes_si = meshsi.nodes;
connect_edge_si = meshsi.connect_edge;

axes(ax);
hold on;

edges_not_hitting_boundary = connect_edge_si(:,3)>=0;

vx = [nodes_si(connect_edge_si(:,1),1),nodes_si(connect_edge_si(:,2),1);nodes_si(connect_edge_si(edges_not_hitting_boundary,2),1),nodes_si(connect_edge_si(edges_not_hitting_boundary,3),1)].';
vy = [nodes_si(connect_edge_si(:,1),2),nodes_si(connect_edge_si(:,2),2);nodes_si(connect_edge_si(edges_not_hitting_boundary,2),2),nodes_si(connect_edge_si(edges_not_hitting_boundary,3),2)].';
h = plot(vx, vy, '-r'); %,nodes_si(:,1),nodes_si(:,2),'.r');
%set(h(1:end-1),'xliminclude','off','yliminclude','off')


if bmarky==1
    % edges with bmark=1
    bm1 = meshsi.bmark_edge(:,1)==1;
    eg1a = meshsi.nodes(meshsi.connect_edge(bm1,1),:);
    eg1b = meshsi.nodes(meshsi.connect_edge(bm1,2),:);   
    eg1 = zeros(size(eg1a,1)*3,2)*NaN; % we need NaN inserted
                                       % between edges to stop
                                       % wrong edges being drawn
    eg1(1:3:end,:) = eg1a;
    eg1(2:3:end,:) = eg1b ;   
    plot(eg1(:,1),eg1(:,2), 'k')

    bm1 = meshsi.bmark_edge(:,2)==1;
    eg1a = meshsi.nodes(meshsi.connect_edge(bm1,2),:);
    eg1b = meshsi.nodes(meshsi.connect_edge(bm1,3),:);   
    eg1 = zeros(size(eg1a,1)*3,2)*NaN; % we need NaN inserted
                                       % between edges to stop
                                       % wrong edges being drawn
    eg1(1:3:end,:) = eg1a;
    eg1(2:3:end,:) = eg1b ;   
    plot(eg1(:,1),eg1(:,2), 'k')
    
    % bmark == 2
    bm2 = meshsi.bmark_edge(:,1)==2;
    eg1a = meshsi.nodes(meshsi.connect_edge(bm2,1),:);
    eg1b = meshsi.nodes(meshsi.connect_edge(bm2,2),:);   
    eg1 = zeros(size(eg1a,1)*3,2)*NaN;
    eg1(1:3:end,:) = eg1a;
    eg1(2:3:end,:) = eg1b ;   
    plot(eg1(:,1),eg1(:,2), 'g')

    bm2 = meshsi.bmark_edge(:,2)==2;
    eg1a = meshsi.nodes(meshsi.connect_edge(bm2,2),:);
    eg1b = meshsi.nodes(meshsi.connect_edge(bm2,3),:);   
    eg1 = zeros(size(eg1a,1)*3,2)*NaN;
    eg1(1:3:end,:) = eg1a;
    eg1(2:3:end,:) = eg1b ;   
    plot(eg1(:,1),eg1(:,2), 'g')
end


hold off;
axis equal;