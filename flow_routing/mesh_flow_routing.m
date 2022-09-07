function [path_lengths, path_nodes] = mesh_flow_routing(dmesh, z, ii_moulin)
% mesh_flow_routing calculates the steepest-descent path lengths for a
% provided triangular mesh and node elevation z, allowing for moulins and
% inhomogenous boundary conditions

tic;

dphic_ds = zeros(dmesh.tri.n_edges, 1);
for ii=1:dmesh.tri.n_edges
    neigh_nodes = dmesh.tri.connect_edge(ii, :);
    z1 = dmesh.tri.nodes(neigh_nodes(1));
    z2 = dmesh.tri.nodes(neigh_nodes(2));
    dphic_ds(ii) = ( (z2 - z1)/dmesh.tri.edge_length(ii));
end

path_lengths = zeros(dmesh.tri.n_elements,1);

for ii=1:dmesh.tri.n_elements
   nodal_z = z(dmesh.tri.connect(ii,:));
   
   [~,el_kk] = min(nodal_z);
   kk_current = dmesh.tri.connect(ii,el_kk);
   path_length = 0;
   path_node = [kk_current];
   while dmesh.tri.bmark(kk_current) <= 0 && ~ismember(kk_current, ii_moulin)
%        zi = z(kk_current);
%        neigh_nodes = dmesh.tri.neigh_node(kk_current,:);
        neigh_edges = dmesh.tri.connect_edge_inv{kk_current};
%        z_neigh = z(neigh_nodes==1);
       gradphi_neigh = dphic_ds(neigh_edges)'.*dmesh.tri.flow_dir{kk_current};
       
       [~,edge_ind] = max(gradphi_neigh); % This is obviously broken...
%        edge_ind
       edge_num = neigh_edges(edge_ind);
       
       all_nodes = dmesh.tri.connect_edge(edge_num, :);
       outlet_ind = all_nodes(all_nodes~=kk_current);
       
       dist = norm(dmesh.tri.nodes(kk_current,:) - dmesh.tri.nodes(outlet_ind,:));
       path_length = path_length + dist;
       kk_current = outlet_ind;
       
       path_node = [path_node, kk_current];
       
   end
   
   path_lengths(ii) = path_length;
   path_nodes{ii} = path_node;
end

toc;
