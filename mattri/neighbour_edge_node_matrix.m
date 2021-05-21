function neigh_edge_node = neighbour_edge_node_matrix(mesh)
% neigh_edge_node = neighbour_edge_node_matrix(mesh)
% 
% Low level function.
%
% Makes a matrix which lists for all edges the neighbouring nodes.
% (n_edges x n_nodes

% This m-file is part of the mattri meshing package.
% Copyright (c) 2014, Mauro A Werder
% All rights reserved.
% Licensed under a BSD 2-Clause License, see LICENCE file

n_edges = mesh.tri.n_edges;
n_nodes = mesh.tri.n_nodes;
neigh_edge_node = sparse(n_edges, n_nodes);
for jj = 1:n_edges
    tnodes = mesh.tri.connect_edge(jj,:);
    neigh_edge_node(jj, tnodes(1)) = 1;
    neigh_edge_node(jj, tnodes(2)) = 1;
end

% $$$ dof = mesh.tri.n_edges;
% $$$ neighlt = sparse(dof, dof);
% $$$ ce = mesh.tri.connect_edge;
% $$$ for jj = 1:dof
% $$$     for kk = 1:dof%jj-1
% $$$         % connect_cell gives how nodes are connected to other nodes via edges
% $$$         if any(ce(jj,:)==kk)
% $$$                 neighlt(jj,kk) = 1;
% $$$         end
% $$$     end
% $$$ end
% $$$ 
% $$$ neigh_edge = neighlt;% + neighlt';
