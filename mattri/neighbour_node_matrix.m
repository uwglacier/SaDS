function neigh_node = neighbour_node_matrix(mesh)
% neigh = neighbour_matrix(mesh)
% 
% Low level function.
% 
% Makes a matrix which lists which nodes neighbour which other nodes.

% This m-file is part of the mattri meshing package.
% Copyright (c) 2014, Mauro A Werder
% All rights reserved.
% Licensed under a BSD 2-Clause License, see LICENCE file

% nodes
dof = mesh.tri.n_nodes;
ce1 = mesh.tri.connect_edge(:,1);
ce2 = mesh.tri.connect_edge(:,2);
neigh_node = sparse(dof, dof);
for jj = 1:dof
    tnodes1 = ce1==jj;
    tnodes2 = ce2==jj;    
    neigh_node(jj, ce2(tnodes1)) = 1;
    neigh_node(jj, ce1(tnodes2)) = 1;
end
