function triangle_write_area(filename, area)
%  triangle_write_area(filename, area)
%
% Low level function
%
% Writes a triangle .area file from an array of triangle maximum areas
%
% filename - name of file (will be overwritten!)
% nodes - x,y coord of nodes
%
% optional:
% comment - comment to be placed in the ele file
% 
% c.f. http://www.cs.cmu.edu/~quake/triangle.area.html

% This m-file is part of the mattri meshing package.
% Copyright (c) 2014, Mauro A Werder
% All rights reserved.
% Licensed under a BSD 2-Clause License, see LICENCE file

fid = fopen(filename, 'w');

% header
n_ele = length(area);

fprintf(fid, '%d\n', [n_ele]);

format_str = strcat('%d %f','\n');
% body
for ii = 1:n_ele
    fprintf(fid, format_str, [ii, area(ii)]);
end

fclose(fid);