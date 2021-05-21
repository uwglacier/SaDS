function triangle_write_meshtri(filename_root, meshtri, comment)
%  triangle_write_meshtri(filename_root, meshtri, comment)
%
% Low level function
%
% From input meshtri it writes triangle files:
% - .node
% - .ele
% - .poly
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
if exist(fileparts(filename_root))~=7
    disp(['Creating output dir: ', fileparts(filename_root)])
    mkdir(fileparts(filename_root))
end

%% elements and nodes are trivial
triangle_write_node([filename_root,'.node'], meshtri.nodes, comment, meshtri.bmark, meshtri.bmark);
triangle_write_ele([filename_root,'.ele'], meshtri.connect, comment);

%% .poly
% prepare

% $$$ indsn = meshtri.bmark>0;
% $$$ indsnn = double(indsn);
% $$$ indsnn(indsn>0) = 1:sum(indsn);
indse = meshtri.bmark_edge>0;
% $$$ 
% $$$ nodes = meshtri.nodes(indsn,:);
% $$$ bmark = meshtri.bmark(indsn,:);

% make the connect edge for just these nodes:
connect_edge = meshtri.connect_edge(indse,:);
%connect_edge = [indsnn(connect_edge(:,1)), indsnn(connect_edge(:,2))];
bmark_edge = meshtri.bmark_edge(indse);

%% write to file


fid = fopen([filename_root,'.poly'], 'w');

% comment
if length(comment)>0
    fprintf(fid, '# %s\n', comment);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first header: for nodes (EMPTY)
fprintf(fid, '%s\n', '# Nodes');
n_nodes = 0;
dim = 2;
n_attr = 0;
n_boundary_markes = 1;

fprintf(fid, '%d %d %d %d\n', [n_nodes, dim, n_attr, n_boundary_markes]);


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