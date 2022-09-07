% optimize_moulins.m attempts to find an optimal distribution of moulins,
% defined by the distribution that maximizes the minimum moulin discharge.
%
% This script works by approximating the equilibrium flux using the simple
% flow routing algorithm from the calc_discharge function. It starts with a
% large number of moulins, placed by randomly selecting interior nodes. It
% then iteratively removes the moulin with the lowest flux, recalculating
% the flow accumulation map to account for the increased flow downstream of
% the removed moulin.

% Add the paths we need
addpath(genpath('../cases/shmip/'))
addpath(genpath('../functions/'))

%% Configuration
meshfile = '../meshes/shmip_refined_mesh_02.mat';
output_file = '../cases/shmip/sensitivity/data/randperm_optimized_moulins_R3.txt';
n_moulins = 50;

%% Calculate gradient of surface elevation along channels
dmesh = load(meshfile);
z_nodes = shmip_elevation(dmesh.tri.nodes);
z_elem = shmip_elevation(dmesh.tri.elements);
melt = shmip_melt(z_elem, 86400*160, false);

dphic_ds = zeros(dmesh.tri.n_edges, 1);
for ii=1:dmesh.tri.n_edges
    neigh_nodes = dmesh.tri.connect_edge(ii, :);
    z1 = dmesh.tri.nodes(neigh_nodes(1));
    z2 = dmesh.tri.nodes(neigh_nodes(2));
    dphic_ds(ii) = ( (z2 - z1)/dmesh.tri.edge_length(ii));
end

%% Generate random moulin locations
z_nodes(dmesh.tri.bmark>0) = Inf;   % No moulins on boundaries
[Y, E] = discretize(z_nodes, 0:250:1500);

moulin_density = [0.6, 0.6, 0.5, 0.4, 0.2, 0.2];
% moulin_counts = [75, 90, 75, 50, 25, 15];
ii_moulin = [];

for binindex=1:length(E)-1
    zindices = find(Y==binindex);
    num_moulin = round(moulin_density(binindex)*length(zindices));
    
    rng(binindex);
%     i_moulin = randi(length(zindices), 1, num_moulin);
    i_moulin = randperm(length(zindices), num_moulin);
    ii_moulin = [ii_moulin; zindices(i_moulin)];
end

%% Optimize position
% By iteratively removing the smallest moulin

% Recalculate the elevation so it doesn't have NaNs (this breaks the
% algorithm)
z_nodes = shmip_elevation(dmesh.tri.nodes);
while length(ii_moulin)>n_moulins
    [~, qm] = calc_discharge(dmesh, z_nodes, dphic_ds, ii_moulin, melt);
    [min_flux, min_index] = min(qm);
    ii_moulin(min_index) = [];
    
    % Print the number of moulins remaining
    disp('Number of moulins remaining:')
    disp(length(ii_moulin))
end

save(output_file, 'ii_moulin', '-ascii')

%% Plot the resulting flux
[qc, qm] = calc_discharge(dmesh, z_nodes, dphic_ds, ii_moulin, melt);
figure
edge_plot(gca, dmesh, qc, palettes('-green-1'), 'none')
