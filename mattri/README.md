mattri: triangular meshes from Matlab/Octave
============================================
These matlab scripts provied an interface to the ['Triangle' meshing
software](http://www.cs.cmu.edu/~quake/triangle.html).

(Aside, if you're in the market for a Matlab mesh generator, you
should have a look at
[DistMesh](http://persson.berkeley.edu/distmesh/).)

Licence
-------

The Matlab interface has a 2-clause BSD licence, see LICENCE.  The
Triangle meshing library (contained in the directory triangle/) is
Copyright by Jonathan Richard Shewchuk, and has a non-commercial
copyleft licence, see triangle/LICENCE.

Installation
============
A Linux & a Windows complied version of Triangle is included in the
sub-directory triangle/, if that does not work, then you need to
compile it yourself: check triangle/README.  Make sure that it is set
to be executable ($chmod u+x triangle).

Demo
====
After compiling Triangle run the demo: `example/demo.m`

One-stop meshing function
=========================
make_mesh.m will make you a mesh from within matlab.  See help of the
function.

Note: running this will create some temporary files in your system
temporary directory in a sub-directory `mesh*`.  Usually these get
deleted but if the program crashes they may linger and have to be
deleted by hand.

Mesh plotting
-------------
Plotting can be done with mesh_plot_tri.m

    >> mesh_plot_tri(gca, mesh)
   
plots both meshes. mesh_plot_tri and mesh_plot_dual only one of them.
mesh_plot_tri.m can also label the nodes and edges.

    >> mesh_plot(gca, mesh)

Mesh datastructure
==================

The generated mesh has the following structure:

    mesh = 
        tri: [1x1 struct]
        [vor: [1x1 struct]]
        [si: [1x1 struct]]

with the triangular mesh sub-structure

    tri = 
               type: 'triangular'		% type
              nodes: [91x2 double]	% coords of nodes
              bmark: [91x1 int32]		% =1 if node is on boundary. Can have more values to make more complex BCs
            connect: [148x3 uint32]	% the connectivity matrix between elements and nodes
       connect_edge: [238x2 uint32]	% the connectivity matrix between edges and nodes
         bmark_edge: [238x1 int32]		% =1 if edge is on boundary
               para: [1x1 struct]		% contains parameters used for the mesh generation
            n_nodes: 127              % number of nodes/vertices
            n_edges: 331              % number of edges
         n_elements: 98			% number of elements (triangles)
               area: [98x1 double]        % area of triangles
        edge_length: [184x1 double]       % length of the edges
     edge_midpoints: [184x2 double]       % midpoint coordinates
    connect_edge_inv: {1x98 cell}          % the inverse of connect_edge
    connect_edge_el: [184x2 int32]	% connectivity matrix between edges and elements
         neigh_node: [127x127 double] % neighbour relations for nodes
    neigh_edge_node: [331x127 double] % neighbour relations for edges

and the dual Voronoi mesh sub-structure

    vor = 
            type: 'Voronoi'		% type
           nodes: [212x2 double]	% coords of nodes
           bmark: [212x1 int32]		% =1 if node is on boundary
    connect_edge: [302x2 uint32]	% the connectivity matrix between edges and nodes
      bmark_edge: [302x1 int32]		% =1 if edge is on boundary
    connect_cell: {91x1 cell}		% connectivity between cells
      bmark_cell: [91x1 int32]		% =1 if cell is on boundary
         connect: {91x1 cell}		% connectivity between cells and edges

NOTE mesh.vor.connect connects cells and edges, not cells and nodes.


Simplex mesh
============
simplex_mesh.m can produce a dual simplex mesh to the triangular
mesh.  Which can be plotted with, you guessed it, mesh_plot_simplex.m

Low level functions
====================
All low level functions are labelled as such in their help.

Geometry specification
----------------------
At the low-level, the geometry of the domain has to be specified in
either a .node or .poly file.  See
[1](https://www.cs.cmu.edu/~quake/triangle.node.html) and
[2](https://www.cs.cmu.edu/~quake/triangle.poly.html). The .poly files can
be used to specify complex geometries (including holes), .node files
work only for convex polygons (as far as I remember).

The two scripts triangle_write_node.m and triangle_write_poly.m can
help with the creation of these files from within matlab.  See their
help.

Mesh generation
---------------
Once a .node or .poly file exists a mesh can be made (but probably
you want to use make_mesh.m). Making a mesh:

    >> mesh = mesh_it('example_mesh/box.node', 0.01, 30)

Where 'input/box.node' needs to be the path to a Triangle '.node' or
'.poly' file; can be found here example_mesh/. The other two
parameters are mesh parameters (max area, and min angle).  Note that
min angle should be <-30 otherwise it will not produce a mesh.


Helper
------
The files triangle_*.m are used to read and write to Triangle specific
files.  In particular if you want to mesh your own geometry run 
triangle_write_poly(output_file_name, coords_of_nodes) 
where coords_of_nodes is a nx2 list of consecutive coordinates.

Mesh reordering
---------------

Triangle returns meshes which are not well ordered which is bad for
computation.  We use an approach from
http://webpages.charter.net/reinerwt/fem.htm: "The one main
neighborhood relation between nodes that was used is that nodes are
neighbored if they belong to the same element. That's certainly
simple. So, a large sparse matrix was put together that contains a one
for each pair of nodes that are in the same element and is zero
otherwise. Using Matlab's subroutine symrcm, a permutation of the
nodes was obtained. This sheme was applied iteratively in the hope to
achieve a slightly better result than with just one iteration."

