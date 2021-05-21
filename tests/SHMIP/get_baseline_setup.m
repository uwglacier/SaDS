function [params, Y0] = get_baseline_setup()
% GET_BASELINE_SETUP sets up the baseline model run for sensitivity tests.
%
% [params, Y0] = get_baseline_setup;
%
% Make any changes to params or Y0, then run as usual:
%
%    run_model(params, Y0, output_path)

% Add paths
addpath(genpath('../../../functions/'))
addpath(genpath('../../../model/'))
addpath(genpath('~/MATLAB Add-Ons/'))
addpath(genpath('../forcing/'))

dmesh = load('data/test_mesh.mat');

params = get_default_params(dmesh);
params.overwrite = true;

%% Setup
% Time-stepping
t = 110*86400;
% tend = t + 2*3600;
tend = 256*86400;
dt = 180;
tt = t:(86400/8):tend;
params.tt = tt;

% params.stats = true;
params.overwrite = true;

params.solver = 'odeRK';
params.solver_opts.dt = dt;

%% Parameters and constants
params.exchange = 'ratio';
params.exchange_ratio = 0.2;
params.channel_model = 'nonlinear';

params.Hmin = 0.1;

params.width = 'ratio';
params.r = 3;   % NEW on 6 April 2021
params.lc = 0;

% Turbulent flow parameterization - much better!
params.alphac = 5/3;
params.betac = 3/2;
params.kc = 10;

params.ks = 1;

% Moulins
ii_moulin = load('data/randperm_optimized_moulins_baseline_02.txt');
params.moulins(ii_moulin) = 1;

params.output_fields = {'m_moulin', 'qc', 'qx_sheet', 'qy_sheet', 'dHcdt', 'dphic_ds', 'dhcdt', 'exchange_frac'};

%% Geometry
z_elements = shmip_elevation(dmesh.tri.elements);
z_nodes = shmip_elevation(dmesh.tri.nodes);
z_edges = shmip_elevation(dmesh.tri.edge_midpoints);

%% Inputs
% Use the SHMIP_MELT function to calculate melt with seasonal + diurnal
% variation
params.ms = @(t) shmip_melt(z_elements, t, true);
params.msc= @(t) shmip_melt(z_edges, t, true);
params.mc = @(t) zeros(dmesh.tri.n_edges, 1);
params.fc = @(t) zeros(dmesh.tri.n_edges, 1);

%% Initial conditions
hc0 = 0.0*ones(dmesh.tri.n_edges,1);
hs0 = 0.0*ones(dmesh.tri.n_elements,1);

% Calculate Hc0
steady = load('outputs/spinup_baseline.mat');
qc = abs(steady.outputs.qc(:, end));
Hc0 = (qc./(params.kc*params.r.*abs(steady.outputs.dphic_ds(:, end).^(params.betac-1)))).^(1/(params.alphac+1));
Hc0(Hc0<0.25) = 0.25;


Y0.hs = hs0;
Y0.zs = z_elements;

Y0.hc = hc0;
Y0.Hc = Hc0;
Y0.phic = z_edges - Hc0;

Y0.Vm = zeros(length(find(params.moulins)),1);

end
