function triangle_write_ele(filename, connect, comment)
%  triangle_write_ele(filename, connect, comment)
%
% Low level function
%
% Writes a triangle .ele file from an array of nodal coordinates
%
% filename - name of file (will be overwritten!)
% nodes - x,y coord of nodes
%
% optional:
% comment - comment to be placed in the ele file
% 
% c.f. http://www.cs.cmu.edu/~quake/triangle.ele.html

% This m-file is part of the mattri meshing package.
% Copyright (c) 2014, Mauro A Werder
% All rights reserved.
% Licensed under a BSD 2-Clause License, see LICENCE file

if nargin < 3
    comment = '';
end

fid = fopen(filename, 'w');

% comment
if length(comment)>0
    fprintf(fid, '# %s\n', comment);
end


% header
n_ele = size(connect,1);
noder_per_tri = 3;
n_attr = 0;

fprintf(fid, '%d %d %d\n', [n_ele, noder_per_tri, n_attr]);

format_str = strcat('%d %d %d %d',repmat(' %f',1,n_attr),'\n');
% body
for ii = 1:n_ele
    fprintf(fid, format_str, [ii, connect(ii,:)]);
end

fclose(fid);