function [path_lengths, path_nodes, qc, qm] = calc_moulin_discharge(dmesh, z, melt, ii_moulin)
% calc_moulin_discharge calculates the steepest-descent path lengths for a
% provided triangular mesh and node elevation z, allowing for moulins and
% inhomogenous boundary conditions

tic;

dphic_ds = zeros(dmesh.tri.n_edges, 1);
for ii=1:dmesh.tri.n_edges
    neigh_nodes = dmesh.tri.connect_edge(ii, :);
%     z1 = dmesh.tri.nodes(neigh_nodes(1));
%     z2 = dmesh.tri.nodes(neigh_nodes(2));
    z1 = z(neigh_nodes(1));
    z2 = z(neigh_nodes(2));
    dphic_ds(ii) = ( (z2 - z1)/dmesh.tri.edge_length(ii));
end

path_lengths = zeros(dmesh.tri.n_elements,1);

qc = zeros(dmesh.tri.n_edges, 1);
qm = zeros(length(ii_moulin), 1);

for ii=1:dmesh.tri.n_elements
   nodal_z = z(dmesh.tri.connect(ii,:));
   
   [~,el_kk] = min(nodal_z);
   kk_current = dmesh.tri.connect(ii,el_kk);
   path_length = 0;
   path_node = [kk_current];
   path_edges = [];
%    ii
   while dmesh.tri.bmark(kk_current) <= 0 && ~ismember(kk_current, ii_moulin)
%        zi = z(kk_current);
%        neigh_nodes = dmesh.tri.neigh_node(kk_current,:);
        neigh_edges = dmesh.tri.connect_edge_inv{kk_current};
%        z_neigh = z(neigh_nodes==1);
       gradphi_neigh = dphic_ds(neigh_edges)'.*dmesh.tri.flow_dir{kk_current};
%             gradphi_neigh = gradphi_neigh(~ismember(path_node(mm), neigh
%        end
       
       for mm=1:length(gradphi_neigh)
           if ismember(neigh_edges(mm), path_edges)
               gradphi_neigh(mm) = nan;
           end
       end
      
       if all(isnan(gradphi_neigh))
           break
       end
       
       [~,edge_ind] = max(gradphi_neigh);
       edge_num = neigh_edges(edge_ind);
       
       all_nodes = dmesh.tri.connect_edge(edge_num, :);
       outlet_ind = all_nodes(all_nodes~=kk_current);
       
       qc(edge_num) = qc(edge_num) + melt(ii)*dmesh.tri.area(ii);
       
       dist = norm(dmesh.tri.nodes(kk_current,:) - dmesh.tri.nodes(outlet_ind,:));
       path_length = path_length + dist;
       kk_current = outlet_ind;
       
       path_node = [path_node, kk_current];
       
       % Find the edge we used
       path_edges = [path_edges, edge_num];
       
   end
   
   path_lengths(ii) = path_length;
   path_nodes{ii} = path_node;
   
   for nn=1:length(ii_moulin)
       if kk_current==ii_moulin(nn)
           qm(nn) = qm(nn) + melt(ii)*dmesh.tri.area(ii);
       end
   end
end

toc;
