function vprime=rhs_sheet_and_channel(t,v,dmesh,params,mode)
% function rhs_sheet_and_channel computes all time derivatives for the
% supraglacial hydrology model. This function is excessively long, but
% this is done on purpose to avoid unnecessarily complicated directory
% structures and dependencies.
%
% vprime = rhs_sheet_and_channel(t, v, dmesh, params, mode)
%
% It's not recommended that you use this function directly. Rather, use the run_model function to control your model runs.
%
% See run_model

% Unpack the state vector
[h_s,z_s,h_c,H_c,phic]=unpack_state_vector(dmesh,v);

% Main variables:
phis = z_s + h_s;                   % Sheet potential scaled by rhow*g

% Channel width is a ratio of channel incision depth
wc = params.r*H_c;

if length(find(isnan(v)))>1e3 % Arbitrary threshold
    if strcmp(mode, 'solver')
        vprime = nan(size(v));
    else
        vprime.phis_x = nan(dmesh.tri.n_elements, 1);
        vprime.phis_y = nan(dmesh.tri.n_elements, 1);

        vprime.dhsdt = nan(dmesh.tri.n_elements, 1);
        vprime.dhcdt = nan(dmesh.tri.n_edges, 1);
        vprime.dHcdt = nan(dmesh.tri.n_edges, 1);
        vprime.dhcdt_area = nan(dmesh.tri.n_edges, 1);

        vprime.dhsdt_div = nan(dmesh.tri.n_elements, 1);
        vprime.dhcdt_div = nan(dmesh.tri.n_edges, 1);

        vprime.qx_edge = nan(dmesh.tri.n_edges, 1);
        vprime.qy_edge = nan(dmesh.tri.n_edges, 1);

        vprime.qx_sheet = nan(dmesh.tri.n_elements, 1);
        vprime.qy_sheet = nan(dmesh.tri.n_elements, 1);

        vprime.qc = nan(dmesh.tri.n_edges, 1);

        vprime.dphic_ds = nan(dmesh.tri.n_edges, 1);
        vprime.phic_node = nan(dmesh.tri.n_nodes, 1);

        vprime.Xi_c = nan(dmesh.tri.n_edges, 1);

        vprime.exchange_frac = nan(dmesh.tri.n_edges, 1);
        vprime.m_moulin = nan(length(find(params.moulins==1)), 1);
    end
else
allow_diss = ones(dmesh.tri.n_edges, 1);
% =========================================================================
%% Sheet model
% Compute the gradient of phi using Green-Gauss Least-Squares approach.
phis_x = zeros(size(h_s));    % x-component of gradient of sheet potential
phis_y = zeros(size(h_s));    % y-component of gradient of sheet potential
for ii=1:dmesh.tri.n_elements
    neigh_els = dmesh.tri.node_stencil_extended{ii};

    % To speed up the model we precalculate these matrices and load them as needed
    A = dmesh.tri.sheet_gradient_matrices{ii};
    W = diag(dmesh.tri.sheet_gradient_weights{ii});

    b = phis(neigh_els) - phis(ii);
    % Solve least squares and store gradient values
    gradvec = (W*A)\(W*b);
    phis_x(ii) = gradvec(1);
    phis_y(ii) = gradvec(2);
end

% Use the potential gradient we just computed to calculate the flux at each
% element centroid
qx_sheet=-params.ks.*abs(h_s).^params.alphas.*sign(h_s).*...
            abs(phis_x).^(params.betas-1).*sign(phis_x);

qy_sheet=-params.ks.*abs(h_s).^params.alphas.*sign(h_s).*...
            abs(phis_y).^(params.betas-1).*sign(phis_y);

% Flux defined on edges (size n_edges x 1) is computed as the mean of
% neighbouring elements, or as the upwind flux
qx_edge = zeros(dmesh.tri.n_edges,1);
qy_edge = zeros(dmesh.tri.n_edges,1);
for ii=1:dmesh.tri.n_elements
    for kk=1:3
        iEdge = dmesh.tri.connect_el_edge(ii, kk);

        % Check if we have already calculated this edge value. This is a
        % simple way to avoid double calculating each edge value. We could
        % alternatively write the look over dmesh.tri.n_edges
        if qx_edge(iEdge)==0 && qy_edge(iEdge)==0
            bmark=dmesh.tri.bmark_el(ii,kk);

            q_x1=qx_sheet(ii);
            q_y1=qy_sheet(ii);
            qn1 = dmesh.tri.nx(ii,kk)*q_x1 + dmesh.tri.ny(ii,kk)*q_y1;

            adj_i=dmesh.tri.connect_el_el(ii,kk);

            if bmark==0 % Interior element
                q_x2=qx_sheet(adj_i);
                q_y2=qy_sheet(adj_i);
                qn2 = -dmesh.tri.nx(ii,kk)*q_x2 - dmesh.tri.ny(ii,kk)*q_y2;

                % Upwind flux method
                if qn1>0 && qn2<0
                    % Element 1 is upwind and 2 is downwind
                    qx = q_x1;
                    qy = q_y1;
                elseif qn1<0 && qn2>0
                    % Element 2 is upwind and 1 is downwind
                    qx = q_x2;
                    qy = q_y2;
                elseif qn1<0 && qn2<0
                    qx = 0;
                    qy = 0;
                else
                    qx = 0.5*(q_x1 + q_x2);
                    qy = 0.5*(q_y1 + q_y2);
                end

            else % Apply boundary conditions
                if bmark==1 % Dirichlet boundary
                    q_x2 = -q_x1;
                    q_y2 = -q_y1;

                elseif bmark==2 % Neumann boundary
                    q_x2 = -q_x1 + 2*params.q_neumann(ii)/dmesh.tri.nx(ii,kk);
                    q_y2 = -q_y1 + 2*params.q_neumann(ii)/dmesh.tri.ny(ii,kk);

                    if isinf(q_x2)
                        q_x2=-q_x1;
                    end

                    if isinf(q_y2)
                        q_y2=-q_y1;
                    end

                elseif bmark==3 % Free-flux boundary

                    % NEW 2021-02-26 TH
                    if qn1>0
                        q_x2 = q_x1;
                        q_y2 = q_y1;
                    else
                        q_x2 = -q_x1;
                        q_y2 = -q_y1;
                    end

                end
                qx=0.5*(q_x1 + q_x2);
                qy=0.5*(q_y1 + q_y2);
            end
            qx_edge(iEdge) = qx;
            qy_edge(iEdge) = qy;
        end
    end
end


%% Sheet-channel conservation scheme
% Partition sheet flux into edges according to the edge normal projection
% onto the flux vector
if params.exchange
    F_hs = scatteredInterpolant(dmesh.tri.elements(:, 1),...
        dmesh.tri.elements(:, 2), h_s);
    hsc = F_hs(dmesh.tri.edge_midpoints);
    exchange_ratio = params.exchange_ratio;
    exchange_frac = 1 - (h_c - (1-exchange_ratio)*H_c)./(exchange_ratio*H_c + hsc);
    exchange_frac(exchange_frac>1)=1;
    exchange_frac(exchange_frac<0)=0;

    % Bound fraction between 0 and 1
    exchange_frac(exchange_frac<0) = 0;
    exchange_frac(exchange_frac>=1) = 1;

    % No exchange into boundary edges
    exchange_frac(dmesh.tri.bmark_edge>0) = 0;
else
    exchange_frac = zeros(dmesh.tri.n_edges, 1);
end

% Initialize exchange term of the channel model equation
dhcdt_exc = zeros(dmesh.tri.n_edges, 1);

% Initialize divergence term of the sheet model equation
dhsdt_div = zeros(dmesh.tri.n_elements, 1);

qn = zeros(dmesh.tri.n_elements, 3);

q_interior = dmesh.tri.nx.*qx_edge(dmesh.tri.connect_el_edge) + ...
                dmesh.tri.ny.*qy_edge(dmesh.tri.connect_el_edge);
qn(dmesh.tri.bmark_el==0) = q_interior(dmesh.tri.bmark_el==0);
qn(dmesh.tri.bmark_el==3) = q_interior(dmesh.tri.bmark_el==3);

qbndry = repmat(params.q_neumann, 1, 3);

qn(dmesh.tri.bmark_el==1) = 0;
qn(dmesh.tri.bmark_el==2) = qbndry(dmesh.tri.bmark_el==2);

% Now calculate mass loss from the elements
for ii=1:dmesh.tri.n_elements
    for kk=1:3
        iEdge = dmesh.tri.connect_el_edge(ii, kk);
        dhcdt_exc(iEdge) = dhcdt_exc(iEdge) + exchange_frac(iEdge).* max(0, qn(ii,kk))./wc(iEdge);
    end
end


% Calculate mass transfer into downhill element
el_edge_flux = qn;
exchange = (1 - exchange_frac(dmesh.tri.connect_el_edge));

% el_edge_flux is an [n_elements x 3] array that keeps track of mass
% exchange across element edges after subtracting the exchange to the
% edges. el_edge_flux(ii,kk) is the mass flux across the kk-th edge of
% element ii, where >0 represents mass leaving the element and <0
% represents mass entering the element
el_edge_flux(qn<0) = exchange(qn<0).*qn(qn<0);

% Divergence term of sheet water thickness rate of change
dhsdt_div = -sum(el_edge_flux.*dmesh.tri.ds, 2)./dmesh.tri.area;

%% Channel model

flowdirs=dmesh.tri.flow_dir;

% Compute the potential on nodes using weighted least-squares method.
phi_node=zeros(dmesh.tri.n_nodes,1);

for kk=1:dmesh.tri.n_nodes
    neigh_edges=dmesh.tri.connect_edge_inv{kk};

    % As with the sheet model we precalculate A and W matrices. Here,
    % A is empty [] if we weren't able to define a large enough stencil to
    % do the least square solution. This can happen at the corners of the
    % domain for example. In this case we just use the mean of the
    % potential on any neighbours
    A = dmesh.tri.node_interp_matrices{kk};
    b = phic(neigh_edges);


    if isempty(A)
        phi_node(kk) = mean(b);
    else
        W = diag(dmesh.tri.node_interp_weights{kk});
        x=(W*A)\(W*b);
        phi_node(kk)=x(1);

        % TH: 18 May
        if params.moulins(kk)==1 && params.correct_moulin_phi
            bh = h_c(neigh_edges);
            xh = (W*A)\(W*bh);
            phi_node(kk) = phi_node(kk) - xh(1);
        end
    end
end

% Compute spatial derivative of channel potential
dphi_ds = (phi_node(dmesh.tri.connect_edge(:, 2)) -...
            phi_node(dmesh.tri.connect_edge(:, 1)))./dmesh.tri.edge_length;


%% Compute all edge fluxes based on h and dphi_ds
q_edges=-params.kc.*wc.*abs(h_c).^(params.alphac).*sign(h_c).*abs(dphi_ds).^(params.betac-1).*sign(dphi_ds);
% Now flux along boundary edges (it is unstable if we allow this)
q_edges(dmesh.tri.bmark_edge>0) = 0;

%% Compute routing of water between nodes
% flux_arr is a sparse n_node x n_edge array represent water flux into
% nodes from edges. This is the scheme we use to ensure mass conservation
flux_arr=0*dmesh.tri.neigh_edge_node';
n_moulins = length(find(params.moulins));
moulin_flux = sparse(dmesh.tri.n_nodes,1);
for jj=1:dmesh.tri.n_nodes
    neigh_edges=dmesh.tri.connect_edge_inv{jj};
    flowdirs_jj=flowdirs{jj};

    % Check if this is a boundary node
    node_bmark=dmesh.tri.bmark(jj);

    % How we compute the balance of fluxes depends on if this is an
    % interior or a boundary node
    if node_bmark==0 % Interior nodes
        % For interior nodes, sum up the flux into the node, and route this
        % out through the most downhill node

        local_dphi = dphi_ds(neigh_edges);

        % Set the gradient to nan on edges that lead to a dirichlet node
        downstream_nodes = dmesh.tri.connect_edge(neigh_edges,:);
        ds_nodes = downstream_nodes';
        downstream_nodes = ds_nodes(ds_nodes~=jj);
        downstream_bmarks = dmesh.tri.bmark(downstream_nodes);

        local_dphi(downstream_bmarks==2) = nan;

        % Find steepest downhill node (dphids>0 means phi increases towards
        % the node = decreases away from the node)
        [max_gradphi,outlet_edge_ii]=max(local_dphi.*flowdirs{jj}');
        outlet_edge=neigh_edges(outlet_edge_ii);

        % If there exists a downhill neighbouring node, route all flux
        % flowing into the node out through the element edge conncting to
        % the downhill node
        if max_gradphi>0
           for ii=1:length(neigh_edges)
                nn=neigh_edges(ii);
                signed_flux=flowdirs_jj(ii)*q_edges(nn);
                if signed_flux>0
                   flux_arr(jj,nn)=signed_flux; % Inward flux

                   if params.moulins(jj)
                        moulin_flux = moulin_flux + sparse(jj, 1, signed_flux, dmesh.tri.n_nodes, 1);
                   else
                        flux_arr(jj,outlet_edge)=flux_arr(jj,outlet_edge)-signed_flux; % Keep track of outward flux
                   end
                end
           end
        else
            % If node is a local minimum potential we maintain zero flux
            % into the node in the channels and let water build up in the
            % channels until it is not longer a minimum. Wo we just need
            % to set q_edges = 0 so there is no extra heat dissipation
            % Xi_c below) here as this leads to runaway water depth
%             q_edges(neigh_edges) = 0;

            if params.moulins(jj)
                for ii=1:length(neigh_edges)
                    nn=neigh_edges(ii);
                    signed_flux=flowdirs_jj(ii)*q_edges(nn);
                    if signed_flux>0
                       flux_arr(jj,nn)=signed_flux; % Inward flux
                            moulin_flux = moulin_flux + sparse(jj, 1, signed_flux, dmesh.tri.n_nodes, 1);
    %                    else % No outward flux from this node - so mass should
    %                    build up?
    %                         flux_arr(jj,outlet_edge)=flux_arr(jj,outlet_edge)-signed_flux; % Keep track of outward flux
                    end
                end
            else
                allow_diss(neigh_edges) = 0;
            end
        end
    % For boundary nodes, the behaviour depends on the boundary type
    elseif node_bmark==1
        % Dirichlet flux boundary conditions. Enforce zero flux through the
        % boundary.
        flux_arr(jj,:)=0;
    elseif node_bmark==2
       % Neumann boundary conditions. We specify flux through the boundary.
       % The prescribed flux is all routed to the most-downhill element
       % edge, or the least uphill if all elements are uphill. If
       % neumann_flux=0, this condition conserves mass in the channels

       dphids_neighs=dphi_ds(neigh_edges);
       dphids_neighs(dmesh.tri.bmark_edge(neigh_edges)~=0)=nan; % Make sure outlet edge is not a boundary edge
       [max_gradphi,outlet_edge_ii]=nanmax(dphids_neighs.*flowdirs{jj}');
       outlet_edge=neigh_edges(outlet_edge_ii);

       for ii=1:length(neigh_edges)
            nn=neigh_edges(ii);
            signed_flux=flowdirs_jj(ii)*q_edges(nn);
            if signed_flux>0 && dmesh.tri.bmark_edge(nn)==0
               flux_arr(jj,nn)=signed_flux;

               flux_arr(jj,outlet_edge)=flux_arr(jj,outlet_edge)-signed_flux;
            end
       end

       % Then finally add in the Neumann prescribed flux
       flux_arr(jj,outlet_edge)=flux_arr(jj,outlet_edge)+params.qN(jj);

      elseif node_bmark==3
       % Free-flux boundary conditions. Here we compute the positive flux
       % as in the interior node case, and virtually route the water out
       % through the boundary
        q_edges(neigh_edges);
       for ii=1:length(neigh_edges)
                nn=neigh_edges(ii);
           if dmesh.tri.bmark_edge(nn)==0 % Only compute for non-boundary edges
                signed_flux=flowdirs_jj(ii)*q_edges(nn);
                if signed_flux>0
                   flux_arr(jj,nn)=signed_flux;
                   if params.moulins(jj)
                       moulin_flux = moulin_flux + sparse(jj, 1, signed_flux, dmesh.tri.n_nodes, 1);
                   end
                end
           end
       end
    end
end

%% Channel divergence term
% Now that we know the routing of water through the channels this is
% actually pretty easy!

% Flux at the first node for each edge. Negative sign means flux entering through this node
q1 = -diag(flux_arr(dmesh.tri.connect_edge(:, 1), :));
% Flux at the second node for each edge. Positive sign means flux leaving through this node
q2 = diag(flux_arr(dmesh.tri.connect_edge(:, 2), :));

% Calculate flux divergence and divergence term of time derivative
dqds = (q2 - q1)./dmesh.tri.edge_length;
dhcdt_div = -dqds./wc;

%% Channel heat dissipation
Xi_c = params.rhow*params.g*abs(q_edges.*dphi_ds).*allow_diss;
dhcdt_Xi = Xi_c/params.rhow/params.L./wc;
dhcdt_Xi(dmesh.tri.bmark_edge>0) = 0;  % No heat dissipation on boundaries

%% Total time derivatives
% dHcdt_melt = (params.mf - 1) * params.msc(t) * (params.rhow/params.rhoi);
dHcdt_melt = (params.mc(t) - params.msc(t))*(params.rhow/params.rhoi);
% dHcdt_melt(H_c <= params.Hmin) = 0;
%     subset = H_c + params.solver_opts.dt*dHcdt_melt < params.Hmin;
%     dHcdt_melt(subset) = (params.Hmin - H_c(subset))./params.solver_opts.dt;

dHc_Xi = (params.rhow/params.rhoi)*dhcdt_Xi;

% Total incision is sum of meltout and viscous heat melt of stream base
dHcdt = 0.5*dHcdt_melt + 0.5*dHc_Xi;
dHc_resid = (params.Hmin - H_c)/params.solver_opts.dt;

if params.regrow_channels = true
    % The good solution:

    subset = (H_c + dHcdt*params.solver_opts.dt) < params.Hmin;
    dHcdt(subset) = dHc_resid(subset);
else
    % Channels can not regrow
    dHcdt(H_c <= params.Hmin) = 0;
end

dwcdt = params.r*dHcdt;
dzcdt = -params.rhow/params.rhoi*params.mc(t);

%% Total water sheet time derivative
% Sheet heat dissipation
Xi_s = params.rhow*params.g*params.Xi_s*abs(qx_sheet.*phis_x + qy_sheet.*phis_y);

% Sum up individual contributions for total time derivative
dhsdt = dhsdt_div + params.ms(t) + Xi_s/params.rhow/params.L;

%% Total channel time derivatives
% Total channel change is sum of individual components
dhcdt_melt = params.mc(t);
msc = params.msc(t);

dhcdt_melt(H_c < params.Hmin) = msc(H_c < params.Hmin);
dhcdt_area = -abs(h_c./wc).*dwcdt;

dhcdt = dhcdt_div + dhcdt_exc + dhcdt_Xi + dhcdt_melt + dhcdt_area + params.fc(t);
dhcdt(dmesh.tri.bmark_edge>0) = 0; % Maintain zero mass in boundary edges

dzsdt = -params.ms(t) * (params.rhow/params.rhoi) - Xi_s/params.rhoi/params.L;

% Ensure no change in domain boundaries
dhcdt(dmesh.tri.bmark_edge>0)=0;

dphicdt = dzcdt + dhcdt - dHcdt;


dVmdt = zeros(n_moulins, 1);
for vv=1:n_moulins
    moulin_indices = find(params.moulins);
    ii_moulin = moulin_indices(vv);
    dVmdt(vv) = moulin_flux(ii_moulin);
end

if strcmp(mode,'solver')
    if strcmp(params.model,'channel')
        dhsdt=0*dhsdt;
        dzsdt=0*dzsdt;
    elseif strcmp(params.model,'sheet')
        dHcdt=0*dHcdt;
        dphicdt=0*dphicdt;
        dhcdt=0*dhcdt;
    end
    vprime=[dhsdt;dzsdt;dhcdt;dHcdt;dphicdt;dVmdt];

else

    vprime.phis_x = phis_x;
    vprime.phis_y = phis_y;

    vprime.dhsdt = dhsdt;
    vprime.dhcdt = dhcdt;
    vprime.dHcdt = dHcdt;
    vprime.dhcdt_area = dhcdt_area;

    vprime.dhsdt_div = dhsdt_div;
    vprime.dhcdt_div = dhcdt_div;

    vprime.qx_edge = qx_edge;
    vprime.qy_edge = qy_edge;

    vprime.qx_sheet = qx_sheet;
    vprime.qy_sheet = qy_sheet;

    vprime.qc = q_edges;

    vprime.dphic_ds = dphi_ds;
    vprime.phic_node = phi_node;

    vprime.Xi_c = Xi_c;

    vprime.exchange_frac = exchange_frac;

    vprime.m_moulin = moulin_flux(find(params.moulins),:);


end

end
end
