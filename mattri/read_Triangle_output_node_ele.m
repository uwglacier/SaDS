function [nodes, bmark, connect] = read_Triangle_output_node_ele(filename, number)
% Reads the output from running Triangle: only nodes and connect
%
% Low level function
%
% If number is not given 1 is assumed (i.e. read files of form *.1.*)

% This m-file is part of the mattri meshing package.
% Copyright (c) 2014, Mauro A Werder
% All rights reserved.
% Licensed under a BSD 2-Clause License, see LICENCE file

if nargin<2
    number = 1;
end

numberstr = ['.', num2str(number)];

[path,filenameroot,extension] = fileparts(filename);
% TRIANGURLAR MESH output
triroot = strcat(path,'/',filenameroot,numberstr);
[nodes,bmark] = triangle_read_node(strcat(triroot,'.node'));
[connect] = triangle_read_ele(strcat(triroot,'.ele'));
