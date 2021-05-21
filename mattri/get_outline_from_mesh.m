function outline  =  get_outline_from_mesh(meshtri)
%  outline  =  get_outline_from_mesh(meshtri)
%
% Returns (x,y) cooridnates of domain boundary.  (only works for simply connected regions.
%
% Todo: return bmark_edge and bmark at outline as well.

% This m-file is part of the mattri meshing package.
% Copyright (c) 2014, Mauro A Werder
% All rights reserved.
% Licensed under a BSD 2-Clause License, see LICENCE file

bedges = find(meshtri.bmark_edge>0);
bedgesn = (1:length(bedges))';
ceb = meshtri.connect_edge(bedges,:);
ce = meshtri.connect_edge;
no = meshtri.nodes;

outline = [];
outline(1,:) = no(ceb(1,1),:);
cur_edge = 1;
for be=1:length(bedgesn)-1
    next_edge = find(ceb(:,1)==ceb(cur_edge,2));
    if length(next_edge)>1
        error('only singly connected regions work')
    end
    outline(be+1,:) = no(ceb(next_edge,1),:);
    cur_edge = next_edge;
end
outline(end+1,:) = outline(1,:);
