function connecting_elements=connect_el_el(sw_mesh)
% Compute element-element connections
connecting_elements=-1*ones(sw_mesh.tri.n_elements,3);

% for ii=1:mesh.tri.n_edges
%     connecting_elements=mesh.tri.connect_edge_el(ii,:);
%     if isempty(find(

for ii=1:sw_mesh.tri.n_elements
    
    for q=1:3
            node1=sw_mesh.tri.connect(ii,q);
        if q==3
            node2=sw_mesh.tri.connect(ii,1);
        else
            node2=sw_mesh.tri.connect(ii,q+1);
        end
            
        
        for jj=1:sw_mesh.tri.n_elements

            if jj~=ii
                neigh_nodes=sw_mesh.tri.connect(jj,:);
                if ismember(node1,neigh_nodes) && ismember(node2,neigh_nodes)
                    connecting_elements(ii,q)=jj;
                end

            end
            
        end
        
        
    end
    
end
% 
% for ix=1:sw_mesh.tri.n_elements
%    % Find neighbouring cell
%    N=1;
%    for mm=1:sw_mesh.tri.n_edges
%        if sw_mesh.tri.connect_edge_el(mm,1)==ix
%            connecting_elements(ix,N)=sw_mesh.tri.connect_edge_el(mm,2);
%            N=N+1;
%        elseif sw_mesh.tri.connect_edge_el(mm,2)==ix
%            connecting_elements(ix,N)=sw_mesh.tri.connect_edge_el(mm,1);
%            N=N+1;
%        end
%    end
% end