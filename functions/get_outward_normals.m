function [mesh_nx,mesh_ny,tang_dist]=get_outward_normals(sw_mesh)
mesh_nx=zeros(sw_mesh.tri.n_elements,3);
mesh_ny=zeros(sw_mesh.tri.n_elements,3);
% norm_dist=zeros(sw_mesh.tri.n_elements,3);
tang_dist=zeros(sw_mesh.tri.n_elements,3);

for ii=1:sw_mesh.tri.n_elements
% Find the specific normal vector connecting two elements, assuming element
% is not on the boundary (for now)
% ii=100;
% Find connecting elements
connecting_els=sw_mesh.tri.connect_el_el(ii,:);

all_nx=zeros(3,1);
all_ny=zeros(3,1);

for lm=1:3
    jj=connecting_els(lm);
% jj=connecting_els(1);
if jj~=-1
nodes_ii=sw_mesh.tri.connect(ii,:);
nodes_jj=sw_mesh.tri.connect(jj,:);

kk=lm;
node_coords=sw_mesh.tri.nodes([nodes_ii,nodes_ii(1)],:);
tx=node_coords(kk+1,1)-node_coords(kk,1);
ty=node_coords(kk+1,2)-node_coords(kk,2);


nx=ty;
ny=-tx;

nnorm=sqrt(nx.^2+ny.^2);

tang_dist(ii,kk)=sqrt(tx.^2+ty.^2);
% norm_dist(ii,kk)=2*sw_mesh.tri.area(ii)/nnorm;

nx=nx/nnorm;
ny=ny/nnorm;

all_nx(kk)=nx;
all_ny(kk)=ny;

% 
% if kk==1
%     delta=node_coords(3,:)-node_coords(2,:);
% else
%     delta=node_coords(kk,:)-node_coords(kk-1,:);
% end
% norm_dist(ii,kk)=abs(delta(1)*nx + delta(2)*ny);
% norm_dist(ii,kk)=norm(sw_mesh.tri.elements(jj,:)-sw_mesh.tri.elements(ii,:));

% Plot it to check if you're right
% node_coords_ii=sw_mesh.tri.nodes(nodes_ii,:);
% node_coords_jj=sw_mesh.tri.nodes(nodes_jj,:);
% 
% cx_ii=mean(node_coords_ii(:,1));
% cy_ii=mean(node_coords_ii(:,2));
% cx_jj=mean(node_coords_jj(:,1));
% cy_jj=mean(node_coords_jj(:,2));

% plot(cx_ii,cy_ii,'ro')
% plot(cx_jj,cy_jj,[cols(lm),'o'])
% plot([cx_ii,cx_ii+0.1*nx],[cy_ii,cy_ii+0.1*ny],cols(lm))
end
end

if ismember(-1,connecting_els)
    isx=find((all_nx.*all_ny)==0);
    for iii=1:length(isx)
        zi=isx(iii);
    nodes_ii=sw_mesh.tri.connect(ii,:);
    node_coords=sw_mesh.tri.nodes([nodes_ii,nodes_ii(1)],:);
    tx=node_coords(zi+1,1)-node_coords(zi,1);
    ty=node_coords(zi+1,2)-node_coords(zi,2);

    nx=ty;
    ny=-tx;
    
    nnorm=sqrt(tx.^2+ty.^2);
    tang_dist(ii,zi)=nnorm;
    nx=nx/nnorm;
    ny=ny/nnorm;
    
%     
%     if zi==1
%         delta=node_coords(3,:)-node_coords(2,:);
%     else
%         delta=node_coords(zi,:)-node_coords(zi-1,:);
%     end
%     norm_dist(ii,zi)=abs(delta(1)*nx + delta(2)*ny);
% norm_dist(ii,zi)=norm(sw_mesh.tri.elements(jj,:)-sw_mesh.tri.elements(ii,:));

    
%     norm_dist(ii,zi)=2*sw_mesh.tri.area(ii)/nnorm;

    
%     all_nx
%     all_ny
%     nx
    all_nx(zi)=nx;
    all_ny(zi)=ny;
    end
end

mesh_nx(ii,:)=all_nx;
mesh_ny(ii,:)=all_ny;
end