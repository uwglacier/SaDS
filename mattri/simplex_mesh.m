function mesh = simplex_mesh(mesh)
%  mesh = simplex_mesh(mesh)
%
% Creates the dual kinda-simplex mesh from a triangulation: i.e. it
% creates a mesh consisting of all gravity lines from the center of
% gravity to the side middle point, plus at the edge it contains
% the triangle lines.
%
% Input: a mesh structure (having mesh.tri), returns the mesh stucture
% with an added field mesh.si containing the kinda-simplex mesh.
% 
% Returns the structure:
% mesh.si.nodes         % the vertices
% mesh.si.edge_midpoints % the midpoints of the edges
% mesh.si.bmark         % boundary mark of the vertices
% mesh.si.connect_edge  % the connection between vertices and edges
% mesh.si.bmark_edge    % boundary mark of the edges
% mesh.si.connect_cell; % connection between cells and cells
% mesh.si.bmark_cell    % boundary mark of the cells
% mesh.si.connect       % connection between edges and cells
% mesh.si.type = 'Kinda simplex mesh'; 
% mesh.si.norm_vecs     % the normal vectors onto the edges: they point to the right into connect(:,1)
% mesh.si.connect_edge_tri % connection between sub-edges and the triangles in the dual mesh. (Trivial in the interior, not so at the edge)
%
% Notes:
% - egdes are defined by two lines: center of grav. to midpoint to
%   other center of grav. Except on the domain boundary where they
%   go: midpoint to triangular mesh node to midpoint. And except
%   near domain boundary: center of grav. to midpoint to -1
%   (i.e. nowhere).
%   -> thus connect_edge has 3 elements per edge
%
%
% References:
% http://www-sop.inria.fr/asclepios/Publications/Herve.Delingette/ijcv99.pdf
% Mishev_Chen_2007: http://onlinelibrary.wiley.com/doi/10.1002/num.20213/abstract

% This m-file is part of the mattri meshing package.
% Copyright (c) 2014, Mauro A Werder
% All rights reserved.
% Licensed under a BSD 2-Clause License, see LICENCE file

nnodes_t = size(mesh.tri.nodes, 1);
nedge_t = size(mesh.tri.connect_edge,1);
nelem_t = size(mesh.tri.connect,1);
nodes_t = mesh.tri.nodes;
c_edge_t = mesh.tri.connect_edge;
c_t = mesh.tri.connect;

% create centers of gravity
c_grav = 1/3*(nodes_t(c_t(:,1),:) + nodes_t(c_t(:,2),:) + nodes_t(c_t(:,3),:));
% create edge midpoints
midpoints = 0.5*(nodes_t(c_edge_t(:,1),:) + nodes_t(c_edge_t(:,2),:));
% get nodes on the boundary
bnodes_t_ids = find(mesh.tri.bmark>0);
bnodes = nodes_t(bnodes_t_ids,:);

% put them in a nodes list
nodes_s = [c_grav;midpoints;bnodes]; % the first part is dual to the
                               % elements. the second part dual to
                               % the edges
bmark_s = [zeros(size(c_grav,1),1);mesh.tri.bmark_edge;mesh.tri.bmark(bnodes_t_ids)];

% make connection between nodes and edges:
% for every edge_t there are two dual simplex edges
connect_edge_s = zeros(nedge_t+length(bnodes_t_ids),3)-5;
bmark_edge_s = zeros(nedge_t+length(bnodes_t_ids),2);
for et = 1:nedge_t
    el_ids = mesh.tri.connect_edge_el(et,:);
    connect_edge_s(et,:) = [el_ids(1), nelem_t+et, el_ids(2)]; % Note: which elmement lies on which side?
    
    % NOTE: we're still missing the edges along the boudary, these
    % will be produced in the next loop.
end

% make connection between edges and cells
% $$$ %%%this was for the connect_s connecting edges to cells (i.e. find edges belonging to a cell
% $$$ bcell_ind = 1;
% $$$ for cell = 1:nnodes_t
% $$$     % find tri-edges connected to node as their id corresponds to the 
% $$$     % simplex edge id
% $$$     edge_s_ids = find(c_edge_t(:,1)==cell |  c_edge_t(:,2)==cell);
% $$$     connect_s{cell} = edge_s_ids;
% $$$ 
% $$$     % Now, if the node/cell is on the boundary then add another
% $$$     % edge to connect_edge_s and append that edge to connect_s
% $$$     if mesh.tri.bmark(cell)>0
% $$$         % the edge_s ids of the two edges ending at the bounary (which equals also to the edge_t id)
% $$$         ab_edge_s_ids = edge_s_ids(connect_edge_s(edge_s_ids,3)==-1);
% $$$         b_nodes_s_ids = connect_edge_s(ab_edge_s_ids,2);
% $$$         % first make the bounary edge
% $$$         connect_edge_s(nedge_t+bcell_ind,:) = [b_nodes_s_ids(1), nelem_t+nedge_t+bcell_ind, b_nodes_s_ids(2)];
% $$$         % and the two bounary marks the of the (two half) edge:
% $$$         bmark_edge_s(nedge_t+bcell_ind, :) = [mesh.tri.bmark_edge(ab_edge_s_ids(1)), mesh.tri.bmark_edge(ab_edge_s_ids(2))];
% $$$         % update connect_s to include the additional edge
% $$$         connect_s{cell} = [connect_s{cell}; nedge_t+bcell_ind];
% $$$         % increment the index
% $$$         bcell_ind = bcell_ind + 1;
% $$$     end
% $$$ end

connect_s = zeros(size(connect_edge_s,1),2);
% all but the edges lying on the domain boundary have a direct
% correspondence the triangular mesh
connect_s(1:nedge_t,:) = mesh.tri.connect_edge;
connect_s_ind = nedge_t+1;  % index for the connect_s entries which
                            % still need to be created

% for the ones on the domain boundary we need a special treatment
% loop over the cells on the boundary:
for ibcell = 1:length(bnodes_t_ids)
    cell = bnodes_t_ids(ibcell);
    % find tri-edges connected to node as their id corresponds to the 
    % simplex edge id
    edge_s_ids = find(c_edge_t(:,1)==cell |  c_edge_t(:,2)==cell);

    % Now, if the node/cell is on the boundary then add another
    % edge to connect_edge_s and append that edge to connect_s
    if mesh.tri.bmark(cell)>0
        % the edge_s ids of the two edges ending at the bounary (which equals also to the edge_t id)
        ab_edge_s_ids = edge_s_ids(connect_edge_s(edge_s_ids,3)==-1);
        b_nodes_s_ids = connect_edge_s(ab_edge_s_ids,2);
        % first make the bounary edge
        connect_edge_s(nedge_t+ibcell,:) = [b_nodes_s_ids(1), nelem_t+nedge_t+ibcell, b_nodes_s_ids(2)];
        % check wether it has the right orientation (on its left is inside the domain)
        vec_on_b = nodes_s(connect_edge_s(nedge_t+ibcell,2),:) - nodes_s(connect_edge_s(nedge_t+ibcell,1),:);
        vec_inside = nodes_s(connect_edge_s(ab_edge_s_ids(1),2),:) - nodes_s(connect_edge_s(ab_edge_s_ids(1),1),:);
        v1x = vec_on_b(1); v1y = vec_on_b(2);
        v2x = vec_inside(1);  v2y = vec_inside(2);
        if atan2( v1x*v2y-v1y*v2x , v1x*v2x+v1y*v2y ) < 0 % change dir
            connect_edge_s(nedge_t+ibcell,:) = [b_nodes_s_ids(2), nelem_t+nedge_t+ibcell, b_nodes_s_ids(1)];
            % and the two bounary marks the of the (two half) edge:
            bmark_edge_s(nedge_t+ibcell, :) = [mesh.tri.bmark_edge(ab_edge_s_ids(2)), mesh.tri.bmark_edge(ab_edge_s_ids(1))];
        else
            % and the two bounary marks the of the (two half) edge:
            bmark_edge_s(nedge_t+ibcell, :) = [mesh.tri.bmark_edge(ab_edge_s_ids(1)), mesh.tri.bmark_edge(ab_edge_s_ids(2))];
        end

        % update connect_s to include the additional edge
        connect_s(connect_s_ind,:) = [cell, -1];
        % increment the index
        connect_s_ind = connect_s_ind + 1;
    else 
        error_should_only_loop_over_boundary
    end

end



% check that all the boundary edges have been created
if any(any(connect_edge_s==-5))
    error_not_all_boundary_edges_created
end

% make connection between cells and cells
for cell = 1:nnodes_t
    connect_cell_s{cell} = int32([c_edge_t(c_edge_t(:,1)==cell, 2); c_edge_t(c_edge_t(:,2)==cell, 1)].');
end

% calcualte area of cells:
area_s = zeros(size(nodes_t,1),1);
for triel = 1:nelem_t     % loop through all triangles and add 1/3
                          % of their area to the 3 cells 
    area_s(c_t(triel,:)) = area_s(c_t(triel,:)) + 1/3*mesh.tri.area(triel);
end

% calculate the normal vectors
norm_vecs = zeros(2,2,size(connect_s,1));
for iedge = 1:size(connect_s,1)
    if connect_edge_s(iedge,3)>0 % normal case not near boundary
        xs = nodes_s(connect_edge_s(iedge,:),1);
        ys = nodes_s(connect_edge_s(iedge,:),2);    
    else % make non-defined node with (NaN,NaN) coords
        xs = [nodes_s(connect_edge_s(iedge,1:2),1); NaN];
        ys = [nodes_s(connect_edge_s(iedge,1:2),2); NaN];
    end        
    norm_vecs_t = [diff(xs).'; diff(ys).'];
    %not we need the not-normalised vectors! norm_vecs = (norm_vecs*diag(1./sqrt(sum(norm_vecs.^2,1))));
    norm_vecs(:,:,iedge) = [0 1;-1 0] *norm_vecs_t; % rotate such that the normal vectors point out of cell connect(iedge,1). I.e. rotate edge vector negative dir.
                                           % [0 1;-1 0] -> negative rotation
                                           % [0 -1;1 0] -> positive rotation
                                           
end

% make connection between edges and triangles
connect_edge_tri = zeros(size(connect_edge_s,1),2);
ind = [1,3];
for iedge = 1:size(connect_s,1)
    for isubedge = 1:2
        if connect_edge_s(iedge,ind(isubedge))==-1 % nothing to do for non-existent edge
            connect_edge_tri(iedge,isubedge) = -1;
            break
        end
        if bmark_edge_s(iedge,isubedge)==0 % normal case: edge in the interior
            connect_edge_tri(iedge,isubedge) = connect_edge_s(iedge,ind(isubedge)); % triangle id equal to node id
        else % the edge lies on the boundary and thus we need to get the right triangle:
            % this search will be slow! xxx
            connect_edge_tri(iedge,isubedge) = connect_edge_s(connect_edge_s(:,2)==connect_edge_s(iedge,ind(isubedge)), 1);
        end
    end
end


% put into structure
mesh.si.nodes = nodes_s; % the vertices
mesh.si.bmark = int8(bmark_s);
mesh.si.connect_edge = int32(connect_edge_s); % the connection between vertices and edges
mesh.si.bmark_edge = int8(bmark_edge_s);
mesh.si.connect_cell = connect_cell_s.'; % connection between cells and cells (Note: already converted to int32 when makeing the cellarray)
mesh.si.bmark_cell = mesh.tri.bmark;
mesh.si.connect = int32(connect_s); % connection between cells and edges
mesh.si.type = 'Kinda simplex mesh';
mesh.si.edge_midpoints = zeros(size(connect_edge_s,1),2,2)*NaN;
mesh.si.edge_midpoints(:,:,1) = 0.5*(nodes_s(connect_edge_s(:,1),:)+nodes_s(connect_edge_s(:,2),:));
non_minus_one_edges = find(connect_edge_s(:,3)>0);
mesh.si.edge_midpoints(non_minus_one_edges,:,2) = 0.5*(nodes_s(connect_edge_s(non_minus_one_edges,2),:)+nodes_s(connect_edge_s(non_minus_one_edges,3),:));
mesh.si.area = area_s;
mesh.si.norm_vecs = norm_vecs;
mesh.si.connect_edge_tri = int32(connect_edge_tri);

%% number of nodes and such
mesh.si.n_cells = size(mesh.tri.nodes,1);
mesh.si.n_edges = size(mesh.si.connect,1);


% figure
% triplot(mesh.tri.connect, mesh.tri.nodes(:,1), mesh.tri.nodes(:,2));
% hold
% plot(c_grav(:,1),c_grav(:,2),'xr')
% plot(midpoints(:,1),midpoints(:,2),'r.')