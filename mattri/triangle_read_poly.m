function [nodes,mark]=triangle_read_poly(filename)
% Reads a triangle .poly file
%
% Low level function
%
% Returns
% nodes - x,y coord of nodes
% mark - node on boundary if ==1
% 
% c.f. https://www.cs.cmu.edu/~quake/triangle.poly.html

% This m-file is part of the mattri meshing package.
% Copyright (c) 2014, Mauro A Werder
% All rights reserved.
% Licensed under a BSD 2-Clause License, see LICENCE file

'not finished, only needed for more complex geometries'
return

% header lines
[n_nodes,dim,n_attr,n_boundary_markes]=textread(filename, '%d%d%d%d', 1)

% rest
[node_num,node_x,node_y,attr,bmark]=textread(filename, '%d%f%f%f%d', -1,...
                                                'commentstyle', 'shell', 'headerlines',1)

nodes = [node_x,node_y];