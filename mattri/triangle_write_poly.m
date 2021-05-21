function triangle_write_poly(filename, nodes, varargin)
% triangle_write_poly(filename, nodes, varargin)
%
% Low level function
%
% Writes a triangle .poly file from a array of nodal coordinates
%
% filename - name of file (will be overwritten!)
% nodes - x,y coord of nodes (do not close the loop by including one coordinate twice!)
%
% optional (in this order):
% comment - comment to be placed in the node file  (default '')
% bmark - to mark whether node is on boundary      (default -1)
% bmark_edge - to mark whether node is on boundary (default -1)
% random_extra_nodes - add these random extra nodes (default none)
%                      Note that these should all lie inside the domain.
% 
% c.f. https://www.cs.cmu.edu/~quake/triangle.poly.html
%
% Note: does not support geometries with holes (even though .poly
% files do)

% This m-file is part of the mattri meshing package.
% Copyright (c) 2014, Mauro A Werder
% All rights reserved.
% Licensed under a BSD 2-Clause License, see LICENCE file

if nargin > 2
    comment = varargin{1};
else
    comment = '';
end
if nargin > 3
    bmark = varargin{2};
else
    bmark = -1;
end
if nargin > 4
    bmark_edge = varargin{3};
else
    bmark_edge = -1;
end
if nargin > 5
    random_extra_nodes = varargin{4};
else
    random_extra_nodes = [];
end

if exist(fileparts(filename))~=7
    disp(['Creating output dir: ', fileparts(filename)])
    mkdir(fileparts(filename))
end


fid = fopen(filename, 'w');

% comment
if length(comment)>0
    fprintf(fid, '# %s\n', comment);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first header: for nodes
fprintf(fid, '%s\n', '# Nodes');
n_nodes = size(nodes,1);
if ~isempty(random_extra_nodes)
    n_nodes_tot = n_nodes+size(random_extra_nodes,1);
else
    n_nodes_tot = n_nodes;
end
dim = 2;
n_attr = 0;
if bmark==-1
    n_boundary_markes = 0;
else
    n_boundary_markes = 1;
end
fprintf(fid, '%d %d %d %d\n', [n_nodes_tot, dim, n_attr, n_boundary_markes]);

format_str = strcat('%d %f %f',repmat(' %f',1,n_attr),repmat(' %d',1,n_boundary_markes),'\n');
% body
if bmark==-1
    for ii = 1:n_nodes
        fprintf(fid, format_str, [ii, nodes(ii,:)]);
    end
else
    for ii = 1:n_nodes
        fprintf(fid, format_str, [double(ii), nodes(ii,:), double(bmark(ii))]);
    end
end
% extra random points
if ~isempty(random_extra_nodes)
    for ii = n_nodes+1:n_nodes+size(random_extra_nodes,1)
        fprintf(fid, format_str, [double(ii), random_extra_nodes(ii-n_nodes,:), 0]);
    end
end    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% second header for segments
n_segments = n_nodes;
if bmark_edge==-1
    n_boundary_markes_seg = 0;
else
    n_boundary_markes_seg = 1;
end
fprintf(fid, '%s\n', '# Segments');
fprintf(fid, '%d %d\n', [n_segments, n_boundary_markes_seg]);

format_str = strcat('%d %d %d',repmat(' %d',1,n_boundary_markes_seg),'\n');
if n_boundary_markes_seg == 0
    for ii = 1:n_segments-1
        fprintf(fid, format_str, [ii, ii, ii+1]);
    end
    fprintf(fid, format_str, [ii+1, ii+1, 1]);
else
    for ii = 1:n_segments-1
        fprintf(fid, format_str, [ii, ii, ii+1,double(bmark_edge(ii))]);
    end
    fprintf(fid, format_str, [ii+1, ii+1, 1,double(bmark_edge(ii+1))]);    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% third header for holes (not implemented)
n_holes = 0;
fprintf(fid, '%s\n', '# Holes');
fprintf(fid, '%d\n', [n_holes]);

for ii = 1:n_holes
    fprintf(fid, '%d %f %f\n', [ii, holes_x(ii), holes_y(ii)]);
end



fclose(fid);