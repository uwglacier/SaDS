function mesh = process_Triangle_output(nodes, bmark, connect, connect_edge, bmark_edge, ...
                                        nodes_vor, bmark_vor, connect_edge_vor, bmark_edge_vor, inf_rays_id, ...
                                        order_mesh_yes, degree_constraint, no_voronoi, trim_vor);
% Low level function

% This m-file is part of the mattri meshing package.
% Copyright (c) 2014, Mauro A Werder
% All rights reserved.
% Licensed under a BSD 2-Clause License, see LICENCE file

if ~exist('no_voronoi', 'var')
    no_voronoi = false;
end
if ~exist('trim_vor', 'var') || isempty(trim_vor)
    trim_vor = true;  % trim the Voronoi to the boundary
end

%% simple quantities
mesh.tri.type = 'triangular';
mesh.tri.nodes = double(nodes);  % coords of nodes
mesh.tri.bmark = double(bmark);  % boundary mark of nodes
mesh.tri.connect = double(connect);  % nodes belonging to triangle
mesh.tri.connect_edge = double(connect_edge);  % nodes beloning to edge
mesh.tri.bmark_edge = double(bmark_edge);  % boundary mark of edge
% parameters used for mesh generation
if exist('meshsize','var')
    mesh.tri.para.max_tri_size = meshsize;  % max triangle size
else
    mesh.tri.para.max_tri_size = [];
end
mesh.tri.para.degree_constraint = degree_constraint;  % minimum
                                                      % angle of triangle
% $$$ mesh.tri.para.in_filename = filename; % input file name
% $$$ mesh.tri.para.out_filename = [path,'/',triroot]; % output filename

% set the derived quantities
mesh.tri = mesh_tri_derived(mesh.tri);

mesh.x_extent = [min(mesh.tri.nodes(:,1)), max(mesh.tri.nodes(:,1))];
mesh.y_extent = [min(mesh.tri.nodes(:,2)), max(mesh.tri.nodes(:,2))];

if ~no_voronoi
    %% VORONOI DUAL MESH

    % INTERLUDE: in the mesh.tri find the connection between edges
    % and elements. 
    connect_edge_el = connect_edge_vor;
    connect_edge_el(inf_rays_id,2) = -1;
    mesh.tri.connect_edge_el = int32(connect_edge_el);
    % interlude ends
    
    if ~all(bmark_edge_vor==0)
        error_something_went_wrong;
    end
    if ~all(bmark_vor==0)
        error_something_went_wrong;
    end
    
    % Voronoi topology, c.f. https://www.cs.cmu.edu/~quake/triangle.topo.html
    % note: Voronoi cell number = node number of trimesh
    connect_cell_vor = cell(size(nodes,1),1);
    connect_vor = cell(size(nodes,1),1);
    % connections of cells:  connect_cell_vor
    % and connections between cells and edges: connect_vor
    for edge_tri_id = 1:size(connect_edge,1) % corresponds to  voronoi edge id
        node1 = connect_edge(edge_tri_id,1);  % corresponds to voronoi cell id
        node2 = connect_edge(edge_tri_id,2);   
        connect_cell_vor{node1}(end+1) = int32(node2); % connection between cells
        connect_cell_vor{node2}(end+1) = int32(node1);
        connect_vor{node1}(end+1) = int32(edge_tri_id); % connection between cells and their edges
        connect_vor{node2}(end+1) = int32(edge_tri_id);    
    end
    % cells at the boundary
    bmark_cell_vor = bmark;


    % VORONOI MESH CUTTING
    % now the Voronoi mesh needs to be processed: 
    % trim it to the edge of the domain and insert necessary nodes and
    % edges.
    if trim_vor
        [nodes_vor, bmark_vor, connect_vor, connect_edge_vor, bmark_edge_vor] = ...
            trim_voronoi(nodes_vor, bmark_vor, connect_vor, bmark_cell_vor, ...
                         connect_cell_vor, connect_edge_vor, bmark_edge_vor, inf_rays_id, mesh.tri);
    else
        [nodes_vor, bmark_vor, connect_edge_vor] = ...
            no_trim_voronoi(nodes_vor, bmark_vor, connect_vor, bmark_cell_vor, connect_cell_vor, ...
                            connect_edge_vor, bmark_edge_vor, inf_rays_id, mesh.tri); 
    end
    
    mesh.vor.type = 'Voronoi';
    mesh.vor.nodes = double(nodes_vor); % coordates of nodes
    mesh.vor.bmark = int8(bmark_vor).'; % boundary mark of nodes
    mesh.vor.connect_edge = int32(connect_edge_vor); % nodes belonging to edges
    mesh.vor.bmark_edge = int8(bmark_edge_vor).'; % boundary marker of edges
    mesh.vor.connect_cell = connect_cell_vor; % connectivity of cells (Note: already converted to int32 when makeing the cellarray)
    mesh.vor.bmark_cell = int8(bmark_cell_vor); % bounary marker of cells
    mesh.vor.connect = connect_vor;  % edges bounding a cell (Note: already converted to int32 when makeing the cellarray)
end 

%%%%%%%%%%%%%%%%
% order the mesh
 if order_mesh_yes
    mesh = order_mesh(mesh);
else % need neighbour information
    mesh.tri.neigh_node = neighbour_node_matrix(mesh);
    mesh.tri.neigh_edge_node = neighbour_edge_node_matrix(mesh); 
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [nodes_vor, bmark_vor, connect_vor, connect_edge_vor, bmark_edge_vor] = trim_voronoi...
        (nodes_vor, bmark_vor, connect_vor, bmark_cell_vor, connect_cell_vor, ...
                      connect_edge_vor, bmark_edge_vor, inf_rays_id, trimesh)
% The Voronoi diagram contains infinite edges: these need to be
% trimmed to the boundary of our domain and additional edges need
% to be genereated.
%
% This function modifies (by appending) nodes_vor, bmark_vor,
% connect_vor, connect_edge_vor and bmark_edge_vor such that the
% Voronoi mesh is bounded by the boundary.

% go through all cells on the domain boundary
boundary_cells_id = find(bmark_cell_vor>0); 
n_boundary_cells = length(boundary_cells_id);

% we need to make new nodes, egdes, connections and such
% Note all these will be appended to the back of the existing 
% matrices to not destroy the implicit topology info contained
% therein, c.f. https://www.cs.cmu.edu/~quake/triangle.topo.html
for ii = 1:n_boundary_cells
    cell_id = boundary_cells_id(ii);
    % the edges belonging to cell_id
    edges_vor_id = connect_vor{cell_id};
    % triangular edges attached to the cell center
    edges_tri_id = [find(trimesh.connect_edge(:,1)==cell_id); find(trimesh.connect_edge(:,2)==cell_id)].';
    edges_on_boundary_tri_id = edges_tri_id(find(trimesh.bmark_edge(edges_tri_id)>0));
    tmp = trimesh.connect_edge(edges_on_boundary_tri_id,:);
    tmp = unique(tmp(:));
    cell_neighbours_along_boundary_id = setdiff(tmp,cell_id);

    % coord of the cell center    
    pt = trimesh.nodes(cell_id,:);
    % coords of neighbouring cell centers lying along the bounary
    pto = trimesh.nodes(cell_neighbours_along_boundary_id,:);

    rays_id = edges_vor_id(find(inf_rays_id(edges_vor_id)==1));

    all_dists = [];
    all_xyouts = [];
    % calculate all intersections
    for rr = 1:length(rays_id)
        xy = nodes_vor(connect_edge_vor(rays_id(rr),1),:);
        dir = [real(connect_edge_vor(rays_id(rr),2)), ...
               imag(connect_edge_vor(rays_id(rr),2))];
        for pp = 1:size(pto,1)
            % check if line is on boundary of domain
            [xyout,dist] = line_ray_intersect(pt.',pto(pp,:).',xy.',dir.');
            all_dists(rr,pp) = dist;
            all_xyouts(rr,pp,:) = xyout;
        end
    end
    
    inds = []; % the indices of the right intersections
    dists = []; % the distances of the right intersections
    xyouts = []; % the coords of the right intersections
    other_cell_ids = [];
    for rr = 1:length(rays_id)
        inds(rr) = find(min(all_dists(rr,:))==all_dists(rr,:));
        % check if no intersection found
        if isempty(inds(rr))
            mesh.tri = trimesh;
            figure;
            mesh_plot_tri(gca,mesh);
            hold on;                
            tt=[p1;p2;p3];
            plot(tt(:,1),tt(:,2),'k');
            plot(xy(1),xy(2),'.r');
            plot([xy(1),xy(1)+dir(1)],[xy(2) , xy(2)+dir(2)],'r');
            axis equal;
            error_no_intersection_found;
        end
        dists(rr) = all_dists(rr,inds(rr));
        xyouts(rr,:) = squeeze(all_xyouts(rr,inds(rr),:));
        other_cell_ids(rr) = cell_neighbours_along_boundary_id(inds(rr));
        % check if both_rays_intersect_same_line
        if length(other_cell_ids)==2 && other_cell_ids(1)==other_cell_ids(2) 
            mesh.tri = trimesh;                
            figure;
            mesh_plot_tri(gca,trimesh);
            hold on;                
            tt=[p1;p2;p3];
            plot(tt(:,1),tt(:,2),'k');
            plot(xy(1),xy(2),'.r');
            plot([xy(1),xy(1)+dir(1)],[xy(2) , xy(2)+dir(2)],'r');
            xy = nodes_vor(connect_edge_vor(rays_id(jj-1),1),:);
            dir = [real(connect_edge_vor(rays_id(jj-1),2)), ...
                   imag(connect_edge_vor(rays_id(jj-1),2))];
            plot(xy(1),xy(2),'.g');
            plot([xy(1),xy(1)+dir(1)],[xy(2) , xy(2)+dir(2)],'g');
            axis equal;
            error_both_rays_intersect_same_line;
        end
        nodes_vor(end+1,:) = xyouts(rr,:).';  % add the new found node
        new_node_id = size(nodes_vor,1);
        bmark_vor(new_node_id) = 1;          % mark it as a node on the
                                             % boundary
       % change the complex direction number into an index to the
       % new node:
        connect_edge_vor(rays_id(rr),2) = new_node_id;
        % set the inf_rays_id to zero for this ray
        inf_rays_id(rays_id(rr)) = 0;
    end

    % add the cell center as a new node
    nodes_vor(end+1,:) = pt;
    id_cell_center_node = size(nodes_vor,1);
    bmark_vor(id_cell_center_node) = 1;          % mark it as a node on the boundary    

    % now add the two lines along the edge:
    % this is non-trivial as above loop may not be run!
    id_hanging_nodes = find_hanging_nodes(connect_vor{cell_id}, connect_edge_vor);
    if length(id_hanging_nodes)~=2
        error('Something is amiss with hanging nodes.  Cause could be a self intersecting boundary.');
    end
    
    % make the two new edges lying along the boundary
    for kk = 1:2
        connect_edge_vor(end+1,:) = [id_cell_center_node;...
                                     id_hanging_nodes(kk)];
        new_edge_id = size(connect_edge_vor,1);
        % it needs to have the same bmark as the corresponding bmark of the
        % edge of the tri-mesh:
        
        % get the vor edge perpendicular to the boudary:
        edge_perp = edges_vor_id(connect_edge_vor(edges_vor_id,1)==id_hanging_nodes(kk)|connect_edge_vor(edges_vor_id,2)==id_hanging_nodes(kk));
        its_bmark = trimesh.bmark_edge(edge_perp);
        
        bmark_edge_vor(new_edge_id) = its_bmark;
        connect_vor{cell_id}(end+1) = new_edge_id;        
    end

    % a last check
    if bmark_cell_vor(cell_id)==0
        error_cell_bmark;
    end
end % for ii = 1:n_boundary_cells

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function singles = find_hanging_nodes(edges, connect_edge_vor)
% finds the nodes of the edges which are only connected to one node
nodes_of_edges = connect_edge_vor(edges,:);
sorted_edges = sort(nodes_of_edges(:));
inds = diff(sorted_edges)==0;
inds_of_doubles = [inds; 0] + [0; inds];
singles = sorted_edges(~inds_of_doubles);

end

function [nodes_vor, bmark_vor, connect_edge_vor] = ...
    no_trim_voronoi(nodes_vor, bmark_vor, connect_vor, bmark_cell_vor, connect_cell_vor, ...
                      connect_edge_vor, bmark_edge_vor, inf_rays_id, trimesh)
% The Voronoi diagram contains infinite edges: these need to be
% trimmed to the boundary of our domain and additional edges need
% to be genereated.
%
% Now this function doesn't trim them to the boudary but just
% shortens them.
%
% This function modifies (by adding) nodes_vor, bmark_vor,
% connect_vor, connect_edge_vor and bmark_edge_vor such that the
% Voronoi mesh is bounded by the boundary.


inf_rays_ind = find(inf_rays_id==1);
all_dists = [];
all_xyouts = [];
% calculate new nodes
for rr = 1:length(inf_rays_ind)
    xy = nodes_vor(connect_edge_vor(inf_rays_ind(rr),1),:);
    dir = [real(connect_edge_vor(inf_rays_ind(rr),2)), ...
        imag(connect_edge_vor(inf_rays_ind(rr),2))];
    nodes_vor(end+1,:) = xy + 10*dir;
    bmark_vor(end+1) = 1;
    connect_edge_vor(inf_rays_ind(rr),2) = size(nodes_vor,1);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xyo,dist] = line_ray_intersect(p1,p2,xy,dir)
% Find intersections of one ray with many line segments, returns
% NaN if there is no intersection.
% 
% format of p1/2 and xyo is [x1 x2 .. xn; y1 y2 .. yn]
%
% input: 
% - p1/2: coordinates of line endpoints
% - xy: (x,y) coordinates of ray origin
% - dir: direction of ray as vector
%
% output: intesection points and distance from xy (in units of dir)

% make sure sizes is right
if size(p1,1)~=2 || size(p2,1)~=2 || ~all(size(xy)==[2,1]) || ~all(size(dir)==[2,1])
    disp('need row verctors as input');
    return
end

% equation for line
% r = a*r + b, 0<=r<=1
a = p2-p1;
b = p1;
num_lines = size(b,2);

size([a, -dir]);

isMatlab = exist('OCTAVE_VERSION', 'var') == 0;
% to suppress:
% Warning: Matrix is singular to working precision. 
if isMatlab
    warning('off', 'MATLAB:singularMatrix'); 
else
    warning('off', 'Octave:singular-matrix-div');
end
for ii=1:num_lines
    sol(:,ii) = [a(:,ii), -dir]\(xy - b(:,ii));
end
if isMatlab
    warning('on', 'MATLAB:singularMatrix'); 
else
    warning('on', 'Octave:singular-matrix-div');
end

dist = sol(2,:);
xyo = repmat(xy,1,num_lines) + dir*dist;

% check if intersection is on line segment
not_on_segment = (sol(1,:)>=1 |sol(1,:)<=0 | dist<0);
dist(not_on_segment) = NaN;
xyo(:,not_on_segment) = NaN;

end
