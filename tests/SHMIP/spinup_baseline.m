% Spinup script for refined mesh baseline case

% Add paths
addpath(genpath('../../../model/'))
addpath(genpath('../../../functions/'))

dmesh = load('../../../meshes/shmip_refined_mesh_02.mat');
output_path = 'outputs/spinup_baseline.mat';

params = get_default_params(dmesh);
params.overwrite = true;

%% Setup
% Time-stepping
t = 182*86400;      % Peak melt
tend = t + 0.5*86400;
dt = 240;
tt=linspace(t,tend,13);
params.tt = tt;

params.solver = 'odeRK';
params.solver_opts.dt = dt;

%% Parameters and constants
params.exchange = 'ratio';
params.exchange_ratio = 0.4;
params.channel_model = 'nonlinear';
params.Xi_s = 1;    % Allow sheet heat dissipation
params.overwrite = true;

params.Hmin = 0.1;

params.width = 'ratio';
params.r = 5;
params.lc = 0;

% Turbulent flow parameterization - much better!
params.alphac = 5/3;
params.betac = 3/2;
params.kc = 10;

params.ks = 1;


% Output fields
params.output_fields = {'m_moulin', 'qc', 'qx_sheet', 'qy_sheet', 'dHcdt', 'dphic_ds', 'dhcdt'};

%% Geometry
z_elements = shmip_elevation(dmesh.tri.elements);
z_nodes = shmip_elevation(dmesh.tri.nodes);
z_edges = shmip_elevation(dmesh.tri.edge_midpoints);

ii_moulin = load('data/randperm_optimized_moulins_baseline_02.txt');
params.moulins(ii_moulin) = 1;

%% Inputs
% Use the SHMIP_MELT function to calculate melt with seasonal + diurnal
% variation
params.ms = @(t) shmip_melt(z_elements, t, false);
params.msc= @(t) shmip_melt(z_edges, t, false);
params.mc = @(t) zeros(dmesh.tri.n_edges, 1);
params.fc = @(t) zeros(dmesh.tri.n_edges, 1);

%% Initial conditions
hc0 = 0.0*ones(dmesh.tri.n_edges,1);
Hc0 = 0.5*ones(dmesh.tri.n_edges, 1);
hs0 = 0.0*ones(dmesh.tri.n_elements,1);

Y0.hs = hs0;
Y0.zs = z_elements;

Y0.hc = hc0;
Y0.Hc = Hc0;
Y0.phic = z_edges - Hc0;

Y0.Vm = zeros(length(find(params.moulins)),1);

%% solver
outs = run_model(params, Y0, output_path);

%% Quick analysis to test outputs


% % Find where the largest moulins are
% ii_moulin = find(outs.params.moulins);
% [m, I] = sort(outs.outputs.m_moulin(:, end), 'descend');
% I_max = I(1:50);
% big_moulins = ii_moulin(I_max);
% moulins = big_moulins;
% moulins = moulins(:);
% baseline_moulins = moulins;
% save('data/moulins_baseline_02.txt', 'moulins', '-ascii')
% 
% % Now setup the sensitivity test M1
% I2 = I(1:100);
% moulins = ii_moulin(I2);
% moulins = moulins(:);
% save('data/moulins_M1.txt', 'moulins', '-ascii')


figure
hold on
edge_plot(gca, dmesh, abs(outs.outputs.qc(:,end)), palettes('-green-1'), 'none')
plot(dmesh.tri.nodes(ii_moulin, 1), dmesh.tri.nodes(ii_moulin, 2), 'ro')

figure
element_plot(dmesh, outs.outputs.hs(:, end))
cmocean('dense')
colorbar
axis image

figure
plot(outs.outputs.tt/86400, outs.outputs.m_moulin')
