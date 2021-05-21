function [mesh, boundary_inds] = make_mesh(boundary_xy, bmark, bmark_edge, maxarea, triareafn, tmpdir,...
                                              boundary_inds2purge, keep_all_boundary_points, voronoi, randomize)
% [mesh, boundary_inds] = make_mesh(boundary_xy, bmark, bmark_edge, maxarea, triareafn, tmpdir,...
%                                   boundary_inds2purge, keep_all_boundary_points, voronoi, randomize);
%
% This function returns a mesh.  By default, it will discharge boundary points which are
% closer together than what the area constraint demands.  This ensures that there are no
% triangles smaller than maxarea, but may dispose of needed points.  (this behaviour can be turned off)
%
% Input:
% - boundary_xy: xy-coordinates of domain boundary (do not close the loop!)
% - bmark:       the boundary mark of above points (needs to be >0, because ==0 means interior)
% - bmark_edge:  the boundary mark of edge after above point (needs to be >0, because ==0 means interior)
% - maxarea:     maximal area of elements (scalar)
% Opitional input (set to [] if not needed):
% - triareafn: function of xy returning an area constraint of that point. (default []) Is
%               evaluated on trianlge midpoints.
% - tmpdir: directory for temporary files.  If none is given a random directory of form mesh* 
%           will be created in the temporary-directory (and deleted, at least when this function completes).
% - boundary_inds2purge: indices of which boundary_inds not to use.  (Usual workflow:
%                        try without boundary_inds2purge, identify which points make the mesh bad
%                        and pass those in)
% - keep_all_boundary_points: if true keep all points in boundary_xy, except the ones
%                             specified in boundary_inds2purge (default = false)
% - voronoi: if true also make Voronoi dual mesh (default false)
% - randomize: randomizes the mesh by adding few points random points.  (default 0). Uses
%              the value for the seed of the generator.  If an array of length two the second element
%              is the number of random points added (default 20).
%
% Output: 
% - mesh
% - boundary_inds -- indices of boundary_xy which were used for the mesh.
%
% Notes: 
% 1) If a triarea function is given then it will itterate until it finds a mesh which does
%    not change anymore.  For the triareafn to work at least one triangle mid point of the
%    coarsest mesh needs to fall onto each refined area or the refined area borders on the 
%    domain boundary.
% 2) It is difficult to get the boundary marks right!  Check them before usage.  See
%    the code for algorithm used to determine them.  bmark_edge is set corresponding to bmark
%    and defaults to 1 at the boundary.
%    See also: http://www.cs.cmu.edu/~quake/triangle.markers.html

% This m-file is part of the mattri meshing package.
% Copyright (c) 2014, Mauro A Werder
% All rights reserved.
% Licensed under a BSD 2-Clause License, see LICENCE file

plot_yes = 0;  % for debugging

degree_constraint = 30;

if ~exist('triareafn', 'var') || isempty(triareafn)
    triareafn = @(xy) maxarea + 0*xy(:,1)*0;  % just return maxarea
    refine = false;
else 
    refine = true;
end
if ~exist('tmpdir', 'var') || isempty(tmpdir)
    deltmp = true;
    while true
        tmpdir = [tempdir(), 'mesh', num2str(randi(100000))];
        if ~exist('tmpdir', 'dir')
            break
        end
    end
else
    deltmp = false;
end
if ~exist('voronoi', 'var') || isempty(voronoi)
    voronoi = true;
end
if ~exist('randomize', 'var') || isempty(randomize)
    randomize = 0;
else
    if length(randomize)==1
        % approximate number of random points inserted:
        n_rand = 20;
    else
        n_rand = randomize(2);
        randomize = randomize(1);
    end
end

% make tmpdir
if ~exist(tmpdir, 'dir')
    stat = mkdir(tmpdir);
    if ~stat
        error(['could not create directory: ', tmpdir]);
    end
end

if ~exist('boundary_inds2purge', 'var') || isempty(boundary_inds2purge)
    boundary_inds2purge = [];
end

if ~exist('keep_all_boundary_points', 'var') || isempty(keep_all_boundary_points)
    keep_all_boundary_points = false;
end

if ~keep_all_boundary_points
    % path length of domain outline
    path_length = [0; cumsum(sqrt(diff(boundary_xy(:,1)).^2 + diff(boundary_xy(:,2)).^2))];
    median_pl = median(diff(path_length));

    % Now figure out which points of boundary_xy to use in the meshing
    boundary_inds = [1]; % first point
    distfn = @(xy1, xy2) sqrt(sum((xy1-xy2).^2));
    edge_lens = sqrt(2*triareafn(boundary_xy));
    max_edge_len = sqrt(2*maxarea);
    for ind = 2:(length(boundary_xy(:,1))-1)
        edge_len = min(edge_lens(ind), max_edge_len);

        %% first check distance between the points for big gaps
        if distfn(boundary_xy(ind,:),boundary_xy(ind-1,:)) > 10*median_pl
            % big gap between ind and ind-1!  Thus include ind-1.
            % (this is not quite satisfactory...)
            boundary_inds(end+1) = ind-1;
        end
        
        %% bmarks: which points to treat specially
        % Heuristic: if the bmark and the two bmark_edge are not the same, then keep the
        % point.  Otherwise it can be disposed of.
        if bmark_edge(ind)~=bmark(ind) || bmark_edge(ind+1)~=bmark(ind) 
            boundary_inds(end+1) = ind;
            continue;
        end

        %% just a normal point
        dist = distfn(boundary_xy(ind,:), boundary_xy(boundary_inds(end),:));
        if dist>=edge_len
            ind = ind;
            boundary_inds(end+1) = ind;
        end
    end
    boundary_inds(end+1) = length(boundary_xy(:,1)); % last point
    boundary_inds = unique(boundary_inds);  % remove duplicates
else
    boundary_inds = 1:size(boundary_xy,1);
end

% now purge boundary_xy of uneeded points:
boundary_inds(boundary_inds2purge) = [];

if plot_yes
    figure();
    plot(boundary_xy(:,1),boundary_xy(:,2),'.-')
    hold
    plot(boundary_xy(boundary_inds,1),boundary_xy(boundary_inds,2), 'r.-')
    legend('orig', 'for mesh')
    drawnow
end

% modify boundary_xy
boundary_xy = boundary_xy(boundary_inds,:);
bmark = bmark(boundary_inds);
bmark_edge = bmark_edge(boundary_inds);

% add random points
if randomize>0
    minx = min(boundary_xy(:,1)); maxx = max(boundary_xy(:,1));
    miny = min(boundary_xy(:,2)); maxy = max(boundary_xy(:,2));
    % seed random number generator
    rng(randomize); 
    rpts = [];
    
    for i=1:n_rand
        while true
            xy = rand(1,2) .* [maxx-minx, maxy-miny] + [minx, miny];
            if inpolygon(xy(1), xy(2), boundary_xy(:,1), boundary_xy(:,2))
                break
            end
        end
        rpts(end+1,:) = xy;
    end
    dist = dist_between_points(rpts);
    rem = dist<max_edge_len;
    rem = triu(rem-diag(diag(rem)));
    rpts(sum(rem,2)>0,:) = [];
    if size(rpts, 2)<1
        error('Choose another value for randomize')
    end
else
    rpts = [];    
end


%% do the meshing
order_mesh_yes = false;
filn = [tmpdir, '/input'];  % filename root
stopping_crit = 0.05;  % if less than so much change in number of nodes then stop refining

triangle_write_poly([filn, '.poly'], boundary_xy, 'Temporary mesh file, can be deleted.', bmark, bmark_edge, rpts);

%mesh = mesh_it([filn, '.poly'], maxarea, degree_constraint, order_mesh_yes);
cur_filn = [filn, '.poly'];
run_Triangle(cur_filn,maxarea,degree_constraint,order_mesh_yes);
[nodes, bmark, connect] = read_Triangle_output_node_ele(cur_filn);


n_nodes = size(nodes,1);
% refine mesh if triareafn is given
meshnr = 0;
if refine
    while true
        meshnr = meshnr + 1;
        tri_midpoints = 1/3*(  nodes(connect(:,1),:) ...
                             + nodes(connect(:,2),:)...
                             + nodes(connect(:,3),:));
        triarea = triareafn(tri_midpoints);
        cur_filn = [filn,'.', num2str(meshnr)];
        triangle_write_area([cur_filn, '.area'], triarea);
        run_Triangle(cur_filn,[],[],1);

        [nodes, bmark, connect] = read_Triangle_output_node_ele(cur_filn, meshnr + 1);
        
        n_nodes_new = size(nodes,1);
        if abs(n_nodes-n_nodes_new)/n_nodes_new < stopping_crit
            % stop if less than stopping_crit difference
            break
        end
        n_nodes = n_nodes_new;
    end
end

%% process mesh
order_mesh_yes = true;
if voronoi
    [nodes, bmark, connect, connect_edge, bmark_edge, ...
     nodes_vor, bmark_vor, connect_edge_vor, bmark_edge_vor, inf_rays_id] ...
        = read_Triangle_output(filn, meshnr + 1, 2);
    
    mesh = process_Triangle_output(nodes, bmark, connect, connect_edge, bmark_edge, ...
                                   nodes_vor, bmark_vor, connect_edge_vor, bmark_edge_vor, inf_rays_id,...
                                   order_mesh_yes, degree_constraint, ~voronoi);
else
    [nodes, bmark, connect, connect_edge, bmark_edge] = read_Triangle_output(filn, meshnr + 1, 2);
    mesh = process_Triangle_output(nodes, bmark, connect, connect_edge, bmark_edge, ...
                                   [], [], [], [], [], order_mesh_yes, degree_constraint, ~voronoi);
end
    
%% tidy up /tmp/
if deltmp
    delete([tmpdir, '/*']);
    rmdir(tmpdir);
end
