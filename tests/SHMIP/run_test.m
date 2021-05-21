% Baseline driver script

addpath(genpath('../../../functions/'))
addpath(genpath('../../../model/'))
addpath(genpath('~/MATLAB Add-Ons/'))

output_path = './outputs/baseline_01.mat';

[params, Y0] = get_baseline_setup;

%% solver
% run_model(params, Y0, output_path);

%% pickup
Y0.Hc = pickup_Hc(output_path);
pickup_output = './outputs/baseline_02_xis.mat';

run_model(params, Y0, pickup_output);
