function mesh = order_mesh(mesh, plot_yes)
% mesh = order_mesh(mesh, plot_yes)
%
% Low level function.
%
% Triangle returns meshes which are not well ordered which is bad for computation.  I use
% an approach from http://webpages.charter.net/reinerwt/fem.htm: "The one main
% neighborhood relation between nodes that was used is that nodes are neighbored if they
% belong to the same element. That's certainly simple. So, a large sparse matrix was put
% together that contains a one for each pair of nodes that are in the same element and is
% zero otherwise. Using Matlab's subroutine symrcm, a permutation of the nodes was
% obtained. This sheme was applied iteratively in the hope to achieve a slightly better
% result than with just one iteration."
%
% Refs:
% http://www.mathworks.com/products/matlab/examples.html?file=/products/demos/shipping/matlab/sparsity.html
% http://people.sc.fsu.edu/~jburkardt/m_src/rcm/rcm.html
%  --> uses rcm
% https://en.wikipedia.org/wiki/Minimum_degree_algorithm
%
% THINK: this will mess up the ordering of nodes in always the same rotation direction!
% And may break my code somehwere (BCs)?
%  --> should be fine as just the nodes get renamed but their order in connect and connect_edge stays
%      Same goes for ordering the edges: their "names" change but nothing else
%
% TODO:
% - order edges as well.  But this may potentially introduce above problem!

% This m-file is part of the mattri meshing package.
% Copyright (c) 2014, Mauro A Werder
% All rights reserved.
% Licensed under a BSD 2-Clause License, see LICENCE file

if nargin<2
    plot_yes = 0;
end

dof = mesh.tri.n_nodes;

%% sort nodes & edges
% find neighbors 
neigh = neighbour_node_matrix(mesh);


% the RCM algo seems fine
per = symrcm(neigh);
% how could I sort edges too?
% per_edge = symrcm(neigh_edge);

% premute the mesh with per
mesh = permute_mesh(mesh, per);

% save the permutation, just in case
mesh.perm.perm_node = per;
mesh.perm.perm_edge = []; %per_edge;

% inverse permutation http://blogs.mathworks.com/loren/2007/08/21/reversal-of-a-sort/
tmp = 1:mesh.tri.n_nodes;
inv_per(per) = tmp;
mesh.perm.inv_perm_node = inv_per;
% $$$ tmp = 1:mesh.tri.n_edges;
% $$$ inv_per_edge(per_edge) = tmp;
% $$$ mesh.perm.inv_perm_edge = inv_per_edge;
mesh.perm.inv_perm_edge = [];

% save the neighbour info
mesh.tri.neigh_node = neigh(per,per);
mesh.tri.neigh_edge_node = neighbour_edge_node_matrix(mesh); %(per_edge,per_edge);

% redo simplex-mesh if existing
if isfield('si',mesh)
    mesh = simplex_mesh(mesh);
end

if plot_yes
    figure
    hold
    spy(neigh)
    spy(neigh(per,per),'r')
    title(['Nodes. Nonzeros: ', num2str(sum(sum(neigh)))])

    figure
    hold
    spy(neigh_edge_node)
    %    spy(neigh_edge(per_edge,per_edge),'r')
    title(['Edges. Nonzeros: ', num2str(sum(sum(neigh_edge)))])
end
