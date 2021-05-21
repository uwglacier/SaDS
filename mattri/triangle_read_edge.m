function [connect_edge,bmark_edge,inf_rays_id]=triangle_read_edge(filename)
% Reads a triangle .edge file
%
% Low level function
%
% retruns:
% connect_edge - the connection matrix of the edges
% bmark -- boundary marks
% inf_rays_id - set to one where there are infinite rays
% 
% c.f. https://www.cs.cmu.edu/~quake/triangle.edge.html

% This m-file is part of the mattri meshing package.
% Copyright (c) 2014, Mauro A Werder
% All rights reserved.
% Licensed under a BSD 2-Clause License, see LICENCE file

% header line
% note that for higher order elements the nodes_per_ele can be
% greater than 3!

fid = fopen(filename);

% header line
tmp = textscan(fid, '%d%d', 1, 'commentstyle', '#');
[n_edge,n_boundary_markes]=deal(tmp{:});
% note that for higher order elements the nodes_per_ele can be
% greater than 3! which is not implemented here.

% make format string
%format_str = strcat('%d%d%d',repmat('%d',1,n_boundary_markes));
 
% read rest of file
%tmp = textscan(fid, format_str, n_edge, 'commentstyle', '#');
tmp = textscan(fid, '%s', 'commentstyle', '#', 'Delimiter', '\n', 'MultipleDelimsAsOne', 1);
fclose(fid);

% parse it. a bit tricky because of voronoi infinite edges
tmp = tmp{1};
[edge_num,tmp] = strtok(tmp);
[connect1,tmp] = strtok(tmp);
connect1 = str2num(char(connect1{:}));
[connect2,tmp] = strtok(tmp);
connect2 = str2num(char(connect2{:}));

inf_rays_id = zeros(length(connect2),1);
% find all infinite rays
if any(connect2==-1) % has infinite rays
    if n_boundary_markes > 0
        error_not_implemented
    end
    infrays = find(connect2==-1);
    inf_rays_id(infrays) = 1;
    [infray_x,tmp] = strtok(tmp);
    infray_x = str2num(char(infray_x{:}));
    [infray_y,tmp] = strtok(tmp);
    infray_y = str2num(char(infray_y{:}));
    for ii = 1:length(infrays)
        ind = infrays(ii);
        connect2(ind) = complex(infray_x(ii), infray_y(ii));
    end
    bmark_edge = [];
else % just a normal edge file
    bmark_edge = [];
    for ii = 1:n_boundary_markes
        [bm,tmp] = strtok(tmp);
        bmark_edge = [bmark_edge, str2num(char(bm{:}))];
    end
end

connect_edge = [connect1,connect2];
inf_rays_id = logical(inf_rays_id);