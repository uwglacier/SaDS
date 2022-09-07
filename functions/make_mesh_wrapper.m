function dmesh2=make_mesh_wrapper(varargin)
% function make_mesh_wrapper wraps the mattri.make_mesh function to
% add attributes that are necessary for element-centered finite volume
% methods. Adds the following attributes to dmesh.tri:
% * elements: [n_elements x 2] array specifying the centroid of each
%       element face. Stored as [centroid_x, centroid_y].
% * connect_el_el: [n_elements x 3] array specifying connectivity betweeen
%       elements. Order of elements is such that connect_el_el(ii,k) is the
%       neighbour to element ii sharing nodes dmesh.tri.connect(ii,k:k+1).
% * ny, nx: [n_elements x 3] array specifying outward unit normal
%       components for each element. Ordered in the same way as connect_el_el,
%       so that normal vector [nx(ii,k),ny(ii,k)] points outward from element
%       ii to its kth neighbour connect_el_el(ii,k).
% * ds: [n_elements x 3] array specifying the length of the edge connecting
%       element ii and its kth neighbour, using the same ordering as above
% * bmark_el: [n_elements x 3] array specifying the bmark for each edge of
%       an elements, using the same ordering as above

% Make the original mesh
dmesh=make_mesh(varargin{:});
dmesh2=dmesh;

% First compute element-element connections
dmesh2.tri.connect_el_el=connect_el_el(dmesh2);

% Compute element centroid
dmesh2.tri.elements=zeros(dmesh2.tri.n_elements,2);
for ii=1:dmesh2.tri.n_elements
   nodes=dmesh2.tri.connect(ii,:);
   node_coords=dmesh2.tri.nodes(nodes,:);
   
   cx=mean(node_coords(:,1));
   cy=mean(node_coords(:,2));
   
   dmesh2.tri.elements(ii,1)=cx;
   dmesh2.tri.elements(ii,2)=cy;
end

% Compute normal vectors and edge lengths
[nx,ny,ds]=get_outward_normals(dmesh2);
dmesh2.tri.ny=ny;
dmesh2.tri.nx=nx;
dmesh2.tri.ds=ds;

bmark_el=zeros(dmesh.tri.n_elements,3);
connect_el_edge=zeros(dmesh.tri.n_elements,3);
% Find boundary mark for each edge of each element
for ii=1:dmesh.tri.n_elements
   nodes=dmesh.tri.connect(ii,:);
   nodes=[nodes,nodes(1)];
   
   edges=zeros(1,3);
   for jj=1:3
       edgei=find(dmesh.tri.connect_edge(:,1)==nodes(jj) & dmesh.tri.connect_edge(:,2)==nodes(jj+1));
       edgej=find(dmesh.tri.connect_edge(:,2)==nodes(jj) & dmesh.tri.connect_edge(:,1)==nodes(jj+1));
       
       if isempty(edgei)
          edgeix=edgej; 
       else
           edgeix=edgei;
       end
       
       edge_bmark=dmesh.tri.bmark_edge(edgeix);
       
       bmark_el(ii,jj)=edge_bmark;
       
       connect_el_edge(ii,jj)=edgeix;
   end
end

dmesh2.tri.bmark_el=bmark_el;
dmesh2.tri.connect_el_edge=connect_el_edge;

%% Loop to find flow directions
% This should be added to the make_mesh_wrapper function
flowdirs=dmesh.tri.connect_edge_inv; % Initialize with proper shape
for kk=1:dmesh.tri.n_nodes
    edges=dmesh.tri.connect_edge_inv{kk};
    dirskk=ones(1,length(edges));
    for jj=1:length(edges)
        edge_nodes=dmesh.tri.connect_edge(edges(jj),:);
        if edge_nodes(1)==kk
            dirskk(jj)=-1;
        end
    end
    flowdirs(kk)={dirskk};
end
dmesh2.tri.flow_dir=flowdirs;

%% Find node stencils
% Compute element-node connections (stencils for Least-Squares methods)
dmesh2.tri.node_stencil_compact={};
dmesh2.tri.node_stencil_extended={};
for ii=1:dmesh2.tri.n_elements
    el_nodes=dmesh2.tri.connect(ii,:);
    stencil=[];
    for jj=1:dmesh2.tri.n_elements
        if jj~=ii % ignoring self-self connections
            neigh_nodes=dmesh2.tri.connect(jj,:); 
            if any(ismember(el_nodes,neigh_nodes))
                stencil=[stencil,jj];
            end
        end
    end

    dmesh2.tri.node_stencil_extended(ii)={stencil};

    dmesh2.tri.node_stencil_compact(ii)={dmesh2.tri.connect_el_el(ii,dmesh2.tri.connect_el_el(ii,:)>0)};
end

dmesh2.tri.edge_stencil={};
for jj=1:dmesh2.tri.n_edges
    edge_nodes=dmesh2.tri.connect_edge(jj,:);
    stencil=[];
    for kk=1:dmesh2.tri.n_elements
        if kk~=jj
            neigh_nodes=dmesh2.tri.connect(kk,:);
            if any(ismember(edge_nodes,neigh_nodes))
                stencil=[stencil,kk];
            end
        end
    end

    dmesh2.tri.edge_stencil(jj)={stencil};
end