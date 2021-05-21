function refmesh = refine_mesh(mesh, triangle_areas)
%  refmesh = refine_mesh(mesh, triangle_areas)
%
% Refines the mesh according to the area constraints in triangle_areas. 
%
% NOTE: this does not order the mesh.  Do it when finished refining with
% order_mesh(refmesh)!
% 
% triangle_areas -- area of elements, if -1 then don't refine.  If scalar is passed in then this is used for all triangles.
%
% BUGS: - does not handle bmarks of new boundary points well.  At the moment they are just
%         set to 2!!
%
% ToDo:
% - allow for more complex topographies by inluding the poly-file 
% - allow input of quantities which will be interpolated onto mesh:
%   - but what to do with channels?
%   - take care when reordering
% - make better bmark stuff
%
% REF:
% http://www.cs.cmu.edu/~quake/triangle.refine.html

% This m-file is part of the mattri meshing package.
% Copyright (c) 2014, Mauro A Werder
% All rights reserved.
% Licensed under a BSD 2-Clause License, see LICENCE file

% make dir for files
tmpdir = tempdir();
randdir = ['aaaa',num2str(randi(100000,1)),'/'];
dirr = [tmpdir, randdir];
mkdir(dirr);

filn = [dirr,'refmesh'];

if mesh.tri.n_elements~=length(triangle_areas) && length(triangle_areas)~=1
    error('wrong triangle_areas input')
end

% write files
no_comment = '';
triangle_write_meshtri(filn, mesh.tri, no_comment);

% $$$ triangle_write_node([filn,'.node'], mesh.tri.nodes, no_comment, mesh.tri.bmark, mesh.tri.bmark);
% $$$ triangle_write_ele([filn,'.ele'], mesh.tri.connect, no_comment);
% $$$ % make a polygon too!
% $$$ % $$$ outline =  mesh.tri.nodes(mesh.tri.bmark>0,:);
% $$$ % $$$ 
% $$$ % $$$ triangle_write_poly([filn, '.poly'], outline, no_comment, mesh.tri.bmark(IA), mesh.tri.bmark_edge(IA));
% $$$ 
% $$$ triangle_write_poly_from_meshtri([filn, '.poly'], mesh.tri);

if length(triangle_areas)~=1
    triangle_write_area([filn,'.area'], triangle_areas);
end


%% meshing
degree_constraint = 30;
run_Triangle(filn);
[nodes, bmark, connect, connect_edge, bmark_edge, ...
          nodes_vor, bmark_vor, connect_edge_vor, bmark_edge_vor, inf_rays_id] =...
    read_Triangle_output(filn);

%% process mesh
order_mesh_yes = 0; % to fix bmarks
refmesh = process_Triangle_output(nodes, bmark, connect, connect_edge, bmark_edge, ...
                               nodes_vor, bmark_vor, connect_edge_vor, bmark_edge_vor, inf_rays_id,...
                               order_mesh_yes,  degree_constraint);

% tidy up
delete([dirr,'*'])
rmdir(dirr)
