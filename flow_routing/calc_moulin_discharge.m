function [path_lengths, path_nodes, qc, qm] = calc_moulin_discharge(dmesh, z, dphic_ds, melt, ii_moulin)
% calc_moulin_discharge calculates channel flux and moulin discharge
% accounting for inhomogeneous boundary conditions and moulins.
% 
% Calculates steepest-descent path lengths for a provided triangular,
% complete flow paths, and discharge into moulins. The code runs with
% DEM depressions, but flow may be routed uphill out of depressions,
% so it is recommended that DEM sinks be filled.
%
% [path_lengths, path_nodes, qc, qm] = calc_moulin_discharge(dmesh, z, dzds, melt, ii_moulin)

tic;

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
    while dmesh.tri.bmark(kk_current) <= 0 && ~ismember(kk_current, ii_moulin)
        neigh_edges = dmesh.tri.connect_edge_inv{kk_current};
        gradphi_neigh = dphic_ds(neigh_edges)'.*dmesh.tri.flow_dir{kk_current};

        % Set gradient to nan for previously used edges so we don't end
        % up with closed loops
        for mm=1:length(gradphi_neigh)
            if ismember(neigh_edges(mm), path_edges)
                gradphi_neigh(mm) = nan;
            end
        end

        % if we have somehow encountered a closed loop, break out of
        % the loop
        if all(isnan(gradphi_neigh))
           break
        end
        
        % Calculate the steepest downhill edge
        [~,edge_ind] = max(gradphi_neigh);
        edge_num = neigh_edges(edge_ind);

        all_nodes = dmesh.tri.connect_edge(edge_num, :);
        outlet_ind = all_nodes(all_nodes~=kk_current);

        % Update edge flux
        qc(edge_num) = qc(edge_num) + melt(ii)*dmesh.tri.area(ii);
        
        % Calculate new path length
        dist = norm(dmesh.tri.nodes(kk_current,:) - dmesh.tri.nodes(outlet_ind,:));
        path_length = path_length + dist;
        
        kk_current = outlet_ind;

        path_node = [path_node, kk_current];
       
       % Find the edge we used
       path_edges = [path_edges, edge_num];
       
   end
   
   path_lengths(ii) = path_length;
   path_nodes{ii} = path_node;
   
   % Add to moulin discharge if relevant
   for nn=1:length(ii_moulin)
       if kk_current==ii_moulin(nn)
           qm(nn) = qm(nn) + melt(ii)*dmesh.tri.area(ii);
       end
   end
end

toc;
