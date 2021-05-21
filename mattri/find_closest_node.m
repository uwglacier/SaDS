function node = find_closest_node(meshtri, loc)
% node = find_closest_node(meshtri, loc);
%
% finds the node closest to the location loc [x,y]

% This m-file is part of the mattri meshing package.
% Copyright (c) 2014, Mauro A Werder
% All rights reserved.
% Licensed under a BSD 2-Clause License, see LICENCE file


% calc distance
dist_sq = sum((meshtri.nodes - repmat(loc, size(meshtri.nodes,1), 1)).^2, 2);

[mind, node] = min(dist_sq);