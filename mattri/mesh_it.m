function mesh = mesh_it(filename, triarea, degree_constraint, order_mesh_yes, verbosity)
% mesh = mesh_it(filename, triarea, degree_constraint, order_mesh_yes, verbosity)
%
% Lower level function, probably you want to use make_mesh.m instead.
%
% Produces a triangular mesh from a .node or .poly file and returns the connectivity
% matrix.  Note that it reorders the mesh for better FEM properties.
%
%
%
% Input
% filename - full path name of .node input file 
% triarea - maximum allowed area
% degree_constraint - maximum allowed angle\
% order_mesh_yes -- whether mesh should be ordered for better FEM numerics (default true)
% verbosity -- how much output on printed on screen
% 
% Output
%    The generated mesh has the following structure:
%mesh = 
%    tri: [1x1 struct]
%    vor: [1x1 struct]
%with the triangular mesh sub-structure
%tri = 
%             type: 'triangular'		% type
%            nodes: [91x2 double]		% coords of nodes
%            bmark: [91x1 int32]		% =1 if node is on boundary
%          connect: [148x3 uint32]	% the connectivity matrix between elements and nodes
%     connect_edge: [238x2 uint32]	% the connectivity matrix between edges and nodes
%       bmark_edge: [238x1 int32]		% =1 if edge is on boundary
%             para: [1x1 struct]		% contains parameters used for the mesh generation
%             area: [98x1 double]        % area of triangles
%      aread_nodes: [148x1 double]       % area around nodes (= area of dual mesh cells)
%      edge_length: [184x1 double]       % length of the edges
%   edge_midpoints: [184x2 double]       % midpoint coordinates
% connect_edge_inv: {1x98 cell}          % the inverse of connect_edge
%          n_nodes: 86                   % number of nodes
%          n_edges: 184     		% number of edges
%       n_elements: 98			% number of elements (triangles)
%  connect_edge_el: [184x2 int32]	% connectivity matrix between edges and elements
%
%and the dual Voronoi mesh sub-structure
%vor = 
%            type: 'Voronoi'		% type
%           nodes: [212x2 double]	% coords of nodes
%           bmark: [212x1 int32]		% =1 if node is on boundary
%    connect_edge: [302x2 uint32]	% the connectivity matrix between edges and nodes
%      bmark_edge: [302x1 int32]		% =1 if edge is on boundary
%    connect_cell: {91x1 cell}		% connectivity between cells
%      bmark_cell: [91x1 int32]		% =1 if cell is on boundary
%         connect: {91x1 cell}		% connectivity between cells and edges
%
%NOTE mesh.vor.connect connects cells and edges, not cells and nodes.
%
% DATATYPES
% nodes - double
% connect* - int32
% bmark*   - int8 range (-128,127)

% This m-file is part of the mattri meshing package.
% Copyright (c) 2014, Mauro A Werder
% All rights reserved.
% Licensed under a BSD 2-Clause License, see LICENCE file


if ~exist('triarea', 'var') || isempty(triarea)
    triarea = [];
end
    
if ~exist('degree_constraint', 'var') || isempty(degree_constraint)
    degree_constraint = 30;
end

if ~exist('order_mesh_yes', 'var') || isempty(order_mesh_yes)
    order_mesh_yes = true;
end

if ~exist('verbosity', 'var') || isempty(verbosity)
    verbosity = 1;
end

run_Triangle(filename, triarea, degree_constraint, verbosity)

[nodes, bmark, connect, connect_edge, bmark_edge, ...
 nodes_vor, bmark_vor, connect_edge_vor, bmark_edge_vor, inf_rays_id] =...
    read_Triangle_output(filename);

%% process mesh
mesh = process_Triangle_output(nodes, bmark, connect, connect_edge, bmark_edge, ...
                               nodes_vor, bmark_vor, connect_edge_vor, bmark_edge_vor, inf_rays_id,...
                               order_mesh_yes, degree_constraint);


