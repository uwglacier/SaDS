function [nodes,bmark,attrs]=triangle_read_node(filename)
% Reads a triangle .node file
%
% Low level function
% 
% Returns
% nodes - x,y coord of nodes
% bmark - node on boundary if ==1
% 
% c.f. https://www.cs.cmu.edu/~quake/triangle.node.html

% This m-file is part of the mattri meshing package.
% Copyright (c) 2014, Mauro A Werder
% All rights reserved.
% Licensed under a BSD 2-Clause License, see LICENCE file

fid = fopen(filename);

% header line
tmp = textscan(fid, '%d%d%d%d', 1, 'commentstyle', '#');
[n_nodes,dim,n_attr,n_boundary_markes]=deal(tmp{:});

% make format string
format_str = strcat('%d%f%f',repmat('%f',1,n_attr),repmat('%d',1,n_boundary_markes));

% read rest of file
tmp = textscan(fid, format_str, n_nodes, 'commentstyle', '#');
% extract from cell array
node_x = tmp{2};
node_y = tmp{3};
if n_attr>0
    attrs = tmp{3+(1:n_attr)};
else
    attrs = [];
end


if n_boundary_markes>0
    bmark = tmp{3+n_attr+1:end};
else
    bmark = [];
end

nodes = [node_x,node_y];

fclose(fid);