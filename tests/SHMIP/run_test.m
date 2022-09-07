% Baseline driver script

set_paths

[params, Y0] = get_baseline_setup;

% Use previous simulation to set initial channel sizes
pickup_output = './outputs/baseline_pickup.mat';
Y0.Hc = pickup_Hc(pickup_output);

model_output = './outputs/baseline.mat';
run_model(params, Y0, model_output);
