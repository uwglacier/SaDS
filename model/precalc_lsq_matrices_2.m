function dmesh2 = precalc_lsq_matrices_2(dmesh)
    % precalc_lsq_matrices_2. Calculate least squares matrices.
    %
    % precalc_lsq_matrices(dmesh2) calculates the matrix and weights for
    %  the least squares problems for the sheet gradient and channel
    %  potential problems. Adds the fields:
    %       sheet_gradient_matrices
    %       sheet_gradient_weights
    %       node_interp_matrices
    %       node_interp_weights
    % to the structure dmesh2.tri
    %
    % See precalc_lsq_matrices. This function creates 2nd order matrices
    % for A

    dmesh2 = dmesh;
    for ii=1:dmesh2.tri.n_elements
        neigh_els = dmesh2.tri.node_stencil_extended{ii};

        % Add secondary neighbours to make sure we have enough points to
        % fit second-order Taylor polynomial
        if length(neigh_els)<5
            for jj=1:length(neigh_els)
                neigh_els = [neigh_els, dmesh2.tri.node_stencil_extended{neigh_els(jj)}];
            end
            neigh_els = neigh_els(neigh_els~=ii);
            neigh_els = unique(neigh_els);

            % Update the stencil to add secondary neighbours
            dmesh2.tri.node_stencil_extended{ii} = neigh_els;
    %             disp(length(neigh_els))
        end

        % Calculate weighting matrix
        dxs = dmesh2.tri.elements(neigh_els,1) - dmesh2.tri.elements(ii,1);
        dys = dmesh2.tri.elements(neigh_els,2) - dmesh2.tri.elements(ii,2);
        inv_dist = 1./sqrt(dxs.^2 + dys.^2);

        % Assemble matrix A
        A = [dxs, dys, dxs.^2, dxs.*dys, dys.^2];

        dmesh2.tri.sheet_gradient_matrices{ii} = A;
        dmesh2.tri.sheet_gradient_weights{ii} = inv_dist.^2;
    end

    % Precompute arrays for calculate node potential
    dmesh2.tri.node_interp_matrices = {};
    dmesh2.tri.node_interp_weights = {};
    for kk=1:dmesh2.tri.n_nodes
        neigh_edges=dmesh2.tri.connect_edge_inv{kk};
%         if length(neigh_edges)>=5
        nodex=dmesh2.tri.nodes(kk,1);
        nodey=dmesh2.tri.nodes(kk,2);

        dx=dmesh2.tri.edge_midpoints(neigh_edges,1)-nodex;
        dy=dmesh2.tri.edge_midpoints(neigh_edges,2)-nodey;

        A=[ones(size(dx)),dx,dy, 0.5*dx.^2, dx.*dy, 0.5*dy.^2];
        W = diag(1./sqrt(dx.^2 + dy.^2));

        if rank(W*A)<6 && length(neigh_edges)>=6
            dmesh2.tri.node_interp_matrices{kk} = [];
            dmesh2.tri.node_interp_weights{kk} = [];
        else
            dmesh2.tri.node_interp_matrices{kk} = A;
            dmesh2.tri.node_interp_weights{kk} = diag(W);
        end
%         end
        if isempty(dmesh2.tri.node_interp_matrices{kk})
            if length(neigh_edges)>=3
%         elseif length(neigh_edges)>=3
                nodex=dmesh2.tri.nodes(kk,1);
                nodey=dmesh2.tri.nodes(kk,2);

                dx=dmesh2.tri.edge_midpoints(neigh_edges,1)-nodex;
                dy=dmesh2.tri.edge_midpoints(neigh_edges,2)-nodey;

                A=[ones(size(dx)),dx,dy];
                W = diag(1./sqrt(dx.^2 + dy.^2));

                if rank(W*A)<3
                    dmesh2.tri.node_interp_matrices{kk} = [];
                    dmesh2.tri.node_interp_weights{kk} = [];
                else
                    dmesh2.tri.node_interp_matrices{kk} = A;
                    dmesh2.tri.node_interp_weights{kk} = diag(W);
                end
            else
                dmesh2.tri.node_interp_matrices{kk} = [];
                dmesh2.tri.node_interp_weights{kk} = [];
            end
        end
    end

end
