function mesh = disorder_mesh(mesh, plot_yes)
% mesh = disorder_mesh(mesh, plot_yes)
%
% inverse of order_mesh: put mesh back like it was made by Triangle.

% This m-file is part of the mattri meshing package.
% Copyright (c) 2014, Mauro A Werder
% All rights reserved.
% Licensed under a BSD 2-Clause License, see LICENCE file


if nargin<2
    plot_yes = 0;
end

dof = mesh.tri.n_nodes;

if plot_yes
    neigh_node_old = mesh.tri.neigh_node;
    neigh_edge_node_old = mesh.tri.neigh_edge_node;
end

% premute the mesh with inverse permutations
mesh = permute_mesh(mesh, mesh.perm.inv_perm_node, mesh.perm.inv_perm_edge);

% find neighbours matrix
mesh.tri.neigh_node = neighbour_node_matrix(mesh);
mesh.tri.neigh_edge_node = neighbour_edge_node_matrix(mesh);

% clear the permutation info
mesh = rmfield(mesh, 'perm');


% redo simplex-mesh if existing
if isfield('si',mesh)
    mesh = simplex_mesh(mesh);
end

if plot_yes
    figure
    hold
    spy(neigh_node_old,'r')
    spy(mesh.tri.neigh_node)
    title(['Nodes. Nonzeros: ', num2str(sum(sum(mesh.tri.neigh_node)))])

    figure
    hold
    spy(neigh_edge_node_old,'r')
    spy(mesh.tri.neigh_edge_node)
    title(['Edges. Nonzeros: ', num2str(sum(sum(mesh.tri.neigh_edge_node)))])

end


