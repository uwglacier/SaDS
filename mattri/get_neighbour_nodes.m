function n_nodes = get_neighbour_nodes(meshtri, node_list)
%  n_nodes = get_neighbour_nodes(meshtri, node_list)
%
% Function which returns a cell array of all the neighbour nodes of
% the elements of node_list

% This m-file is part of the mattri meshing package.
% Copyright (c) 2014, Mauro A Werder
% All rights reserved.
% Licensed under a BSD 2-Clause License, see LICENCE file


n_nodes = {};
for nn = 1:length(node_list)
    n1 = meshtri.connect_edge(find(meshtri.connect_edge(:,1)==node_list(nn)),2);
    n2 = meshtri.connect_edge(find(meshtri.connect_edge(:,2)==node_list(nn)),1);    
    n_nodes{nn} = [n1;n2];
end
