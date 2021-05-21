function triangle_write_poly_from_meshtri(filename, meshtri, comment)
% triangle_write_poly(filename, meshtri)
%
% Low level function
%
% Writes a triangle .poly file from a the boundary of a meshtri.
%
% filename - name of file (will be overwritten!)
% meshtri  - a triangulation 
%
% c.f. https://www.cs.cmu.edu/~quake/triangle.node.html

% This m-file is part of the mattri meshing package.
% Copyright (c) 2014, Mauro A Werder
% All rights reserved.
% Licensed under a BSD 2-Clause License, see LICENCE file

if ~exist('comment','var')
    comment = '';
end

%% prepare

indsn = meshtri.bmark>0;
indsnn = double(indsn);
indsnn(indsn>0) = 1:sum(indsn);
indse = meshtri.bmark_edge>0;

nodes = meshtri.nodes(indsn,:);
bmark = meshtri.bmark(indsn,:);

% make the connect edge for just these nodes:
connect_edge = meshtri.connect_edge(indse,:);
connect_edge = [indsnn(connect_edge(:,1)), indsnn(connect_edge(:,2))];
bmark_edge = meshtri.bmark_edge(indse);

%% write to file

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
dim = 2;
n_attr = 0;
n_boundary_markes = 1;

fprintf(fid, '%d %d %d %d\n', [n_nodes, dim, n_attr, n_boundary_markes]);

format_str = strcat('%d %f %f',repmat(' %f',1,n_attr),repmat(' %d',1,n_boundary_markes),'\n');
% body
for ii = 1:n_nodes
    fprintf(fid, format_str, [double(ii), nodes(ii,:), double(bmark(ii))]);
end    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% second header for segments
n_segments = size(connect_edge,1);
n_boundary_markes_seg = 1;

fprintf(fid, '%s\n', '# Segments');
fprintf(fid, '%d %d\n', [n_segments, n_boundary_markes_seg]);

format_str = strcat('%d %d %d',repmat(' %d',1,n_boundary_markes_seg),'\n');
for ii = 1:n_segments
    fprintf(fid, format_str, [ii, connect_edge(ii,1), connect_edge(ii,2),double(bmark_edge(ii))]);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% second header for holes (NONE)
n_holes = 0;
fprintf(fid, '%s\n', '# Holes');
fprintf(fid, '%d\n', [n_holes]);

for ii = 1:n_holes
    fprintf(fid, '%d %f %f\n', [ii, holes_x(ii), holes_y(ii)]);
end


fclose(fid);