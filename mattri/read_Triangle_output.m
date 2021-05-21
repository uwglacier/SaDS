function [nodes, bmark, connect, connect_edge, bmark_edge, ...
          nodes_vor, bmark_vor, connect_edge_vor, bmark_edge_vor, inf_rays_id] =...
    read_Triangle_output(filename, number, verbosity)
% Reads the output from running Triangle
%
% Low level function
%
% If number is not given 1 is assumed (i.e. read files of form *.1.*)

% This m-file is part of the mattri meshing package.
% Copyright (c) 2014, Mauro A Werder
% All rights reserved.
% Licensed under a BSD 2-Clause License, see LICENCE file

if ~exist('number', 'var') || isempty(number)
    number = 1;
end

if ~exist('verbosity', 'var') || isempty(verbosity)
    verbosity = 2;
end

numberstr = ['.', num2str(number)];

[path,filenameroot,extension] = fileparts(filename);
% TRIANGURLAR MESH output
triroot = strcat(path,'/',filenameroot,numberstr);
[nodes,bmark] = triangle_read_node(strcat(triroot,'.node'));
[connect] = triangle_read_ele(strcat(triroot,'.ele'));
[connect_edge,bmark_edge] = triangle_read_edge(strcat(triroot,'.edge'));

if nargout>6
    % VORONOI DUAL MESH output
    vorroot = strcat(path,'/',filenameroot,numberstr,'.v');
    [nodes_vor, bmark_vor] = triangle_read_node(strcat(vorroot,'.node'));
    [connect_edge_vor, bmark_edge_vor, inf_rays_id] = triangle_read_edge(strcat(vorroot,'.edge'));
end

if verbosity>1
    disp(['Mesh with ', num2str(size(nodes,1)), ' nodes.'])
    disp(' ')
end
