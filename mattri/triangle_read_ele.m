function [connect]=triangle_read_ele(filename)
% Reads a triangle .ele file
%
% Low level function
%
% Returns
% nodes - x,y coord of nodes
% mark - node on boundary if ==1
% 
% c.f. https://www.cs.cmu.edu/~quake/triangle.ele.html

% This m-file is part of the mattri meshing package.
% Copyright (c) 2014, Mauro A Werder
% All rights reserved.
% Licensed under a BSD 2-Clause License, see LICENCE file

fid = fopen(filename);

% header line
tmp = textscan(fid, '%d%d%d', 1, 'commentstyle', '#');
[n_ele,nodes_per_ele,n_attr]=deal(tmp{:});
% note that for higher order elements the nodes_per_ele can be
% greater than 3! which is not implemented here.
if nodes_per_ele > 3
    'Not working for higher order elements.'
    connect = 0;
    return
end

% make format string
format_str = strcat('%d%d%d%d',repmat('%f',1,n_attr));

% read rest of file
tmp = textscan(fid, format_str, n_ele, 'commentstyle', '#');
% extract from cell array
connect1 = tmp{2};
connect2 = tmp{3};
connect3 = tmp{4};

connect = [connect1,connect2,connect3];
fclose(fid);

