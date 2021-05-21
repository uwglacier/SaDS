function triangle_write_node(filename, nodes, comment, bmark, attrs)
% triangle_write_node(filename, nodes, comment, bmark, attrs)
%
% Low level function
%
% Writes a triangle .node file from a array of nodal coordinates
%
% filename - name of file (will be overwritten!)
% nodes - x,y coord of nodes
%
% optional:
% comment - comment to be placed in the node file
% bmark - to mark whether node is on boundary
% attrs  - other floating point attrs
% 
% c.f. https://www.cs.cmu.edu/~quake/triangle.node.html

% This m-file is part of the mattri meshing package.
% Copyright (c) 2014, Mauro A Werder
% All rights reserved.
% Licensed under a BSD 2-Clause License, see LICENCE file

if ~exist('comment','var') || isempty(comment)
    comment = '';
end

if ~exist('bmark','var') || isempty(bmark)
    bmark = -1;
end

if ~exist('attrs','var') || isempty(attrs)
    attrs = -1;
end


fid = fopen(filename, 'w');

% comment
if length(comment)>0
    fprintf(fid, '# %s\n', comment);
end

% header
n_nodes = size(nodes,1);
dim = 2;
if bmark==-1
    n_boundary_markes = 0;
else
    n_boundary_markes = 1;
end
if attrs==-1
    n_attrs = 0;
else
    n_attrs = size(attrs,2);
end

fprintf(fid, '%d %d %d %d\n', [n_nodes, dim, n_attrs, n_boundary_markes]);

format_str = strcat('%d %f %f',repmat(' %f',1,n_attrs),repmat(' %d',1,n_boundary_markes),'\n');
% body
if bmark==-1 & attrs==-1
    for ii = 1:n_nodes
        fprintf(fid, format_str, [ii, nodes(ii,:)]);
    end
elseif attrs==-1
    for ii = 1:n_nodes
        fprintf(fid, format_str, [double(ii), nodes(ii,:), double(bmark(ii,:))]);
    end
elseif bmark==-1
    for ii = 1:n_nodes
        fprintf(fid, format_str, [double(ii), nodes(ii,:), double(attrs(ii,:))]);
    end
else
    for ii = 1:n_nodes
        fprintf(fid, format_str, [double(ii), nodes(ii,:), double(attrs(ii,:)), double(bmark(ii,:))]);
    end
    
end    

fclose(fid);