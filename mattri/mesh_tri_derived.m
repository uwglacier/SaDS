function meshtri = mesh_tri_derived(meshtri)
% meshtri = mesh_tri_derived(meshtri)
%
% Low level function!
%
% Calculates the derived quantites not directly gotten from the
% Triangle output. This is also useful if eg one dimension is rescaled.
%
% Calculates:
% mesh.tri.area
% mesh.tri.edge_length
% mesh.tri.edge_midpoints
% mesh.tri.n_nodes = size(mesh.tri.nodes,1);
% mesh.tri.n_edges = size(mesh.tri.connect_edge,1);
% mesh.tri.n_elements = size(mesh.tri.connect,1);
% mesh.tri.connect_edge_inv: cell array containing the inverse of connect_edge
% mesh.tri.area_nodes: the area around a node (as in the area of the dual mesh cell)

% This m-file is part of the mattri meshing package.
% Copyright (c) 2014, Mauro A Werder
% All rights reserved.
% Licensed under a BSD 2-Clause License, see LICENCE file

%% number of nodes and such
meshtri.n_nodes = size(meshtri.nodes,1);
meshtri.n_edges = size(meshtri.connect_edge,1);
meshtri.n_elements = size(meshtri.connect,1);

% area of elements http://en.wikipedia.org/wiki/Triangle#Using_coordinates
area = zeros(size(meshtri.connect,1),1);
for el =1:size(meshtri.connect,1)
    area(el) = 1/2*det([meshtri.nodes(meshtri.connect(el,:),:), ones(3,1)]);
end
meshtri.area = area;

% length of edges
meshtri.edge_length = sqrt(sum((meshtri.nodes(meshtri.connect_edge(:,2),:)-meshtri.nodes(meshtri.connect_edge(:,1),:)).^2,2));
meshtri.edge_midpoints = 0.5*(meshtri.nodes(meshtri.connect_edge(:,1),:) + meshtri.nodes(meshtri.connect_edge(:,2),:));

% "inverse" of connect_edge: given a node, which edges connect to it, to front and back.
for nd =1:meshtri.n_nodes
    tmp1 = (sort(find(meshtri.connect_edge(:,1)==nd)))';
    tmp2 = (sort(find(meshtri.connect_edge(:,2)==nd)))';
    connect_edge_inv{nd} = [tmp1, tmp2];
end
meshtri.connect_edge_inv = connect_edge_inv;

%% calcualte area of cells around nodes (cells as for a simplex mesh):
area_nodes = zeros(meshtri.n_nodes,1);
for el = 1:meshtri.n_elements     % loop through all triangles and add 1/3
                          % of their area to the 3 cells 
    area_nodes(meshtri.connect(el,:)) = area_nodes(meshtri.connect(el,:)) + 1/3*meshtri.area(el);
end
meshtri.area_nodes = area_nodes;