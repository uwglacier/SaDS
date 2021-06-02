% Baseline driver script

set_paths

% This sets the parameters
[params, Y0] = get_baseline_setup;
dmesh = params.dmesh;

% Modify the time stepping for the spinup case
t0 = 182*86400;     % Start time
t1 = 182.5*86400;   % End time
dT = 3600;          % Output every hour
params.tt = t0:dT:t1;

% Use forcing with no diurnal variation
z_elements = shmip_elevation(dmesh.tri.elements);
z_edges = shmip_elevation(dmesh.tri.edge_midpoints);

params.ms = @(t) shmip_melt(z_elements, t, false);
params.msc = @(t) shmip_melt(z_edges, t, false);

% Where to save outputs
pickup_output = './outputs/baseline_pickup.mat';

% This line runs the model
run_model(params, Y0, pickup_output);
