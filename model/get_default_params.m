function params = get_default_params(dmesh)
% get_default_params returns a structure with default model parameters.
%
% params = get_default_params(dmesh) gets the default parameter values for
% the mesh dmesh.
%
% NOTE that there are some parameters not set here that must be set by the
% user:
%     params.tt
%     params.solver_opts*
%     params.moulins
%     params.ms
%     params.msc
%     params.mc
%     params.fc
%
% *The fields of params.solver_opts depend on the solver:
%   solver          | Fields
%   --------------- | ----------------------------------------
%   odeRK           | params.solver_opts.dt
%   Built-in (ode45)| params.solver_opts = odeset( ... )
%
% See also validate_params

%% Physical constants -----------------------------------------------------
% You probably don't want to change these!
params.rhow = 1000;     % Density of water (kg/m^3)
params.rhoi = 850;      % Density of near-surface ice (kg/m^3)
params.L = 334e3;       % Latent heat of water (J/kg)
params.g = 9.81;        % Gravity (m/s^2)

%% Tunable parameters -----------------------------------------------------
% These are the "knobs" you can turn to tune the model
params.alphac = 5/3;    % Channel flow exponents
params.betac = 3/2;     % Channel flow exponents
params.kc = 10;         % Channel hydraulic conductivity

params.alphas = 5/4;    % Sheet flow exponents
params.betas = 3/2;     % Sheet flow exponents
params.ks = 1.0;        % Sheet hydraulic conductivity

params.wc = 0.1*ones(dmesh.tri.n_edges, 1);        % Width of channels (m)
params.r = 5;           % Ratio of channel width to incision depth (-)

params.Hmin = 1e-3;     % Minimum incision depth (m) before channel melts at same rate as sheet

params.exchange_ratio = 0.2;    % See params.exchange note below
%% Moulins ----------------------------------------------------------------
% Specify moulins using a sparse array of size (n_nodes, 1). Nodes set to 1
% are moulins
params.moulins = sparse(dmesh.tri.n_nodes, 1);

%% Solver parameters ------------------------------------------------------
% Choose ODE solver from:
% 'odeRK': fourth-order explicit Runge-Kutta method
% 'ode45': built-in matlab ode45 method
% 'odeXX': Any other built-in matlab ode method
params.solver = 'odeRK';
params.solver_opts = struct;    % User must populate this structure

%% Time-stepping
% Note - user will always want to specify this themselves
% Specified as
% params.solver_opts.dt = XX

%% Switches ---------------------------------------------------------------
params.model = 'coupled';       % 'sheet' for sheet only, 'channel' for channel only, 'coupled' for both
params.Xi_s = true;             % Set to 1 (or true) for heat dissipation in sheet
params.overwrite = false;       % Overwrite existing model output files model. 'linear' or 'nonlinear'.
params.correct_moulin_phi = false;

params.stop_on_warning = false; % Switch to raise an error if parameter validation fails

% Control behaviour of small channels. If false, sets
% dHcdt(Hc<=params.Hmin) = 0. If true, allows channels to regrow if
% frictional melt outpaces surface lowering
params.regrow_channels = false;

% params.exchange: Control the mass exchange between elements and edges.
% Possible values are {false, 'linear', 'ratio', 'old'}. false does not
% allow mass exchange between the systems. 'linear' interpolates between 1
% (all mass) when channel is dry to 0 (no mass) when channel is full.
% 'ratio' is similar but interpolates over a ratio of the full channel
% depth.
params.exchange = true; % TH 2021-05-04 changed to true/false

%% Inputs/boundary conditions ---------------------------------------------

% Neumann boundary flux for elements
params.q_neumann = zeros(dmesh.tri.n_elements,1);

% Neumann boundary flux for nodes
params.qN = zeros(dmesh.tri.n_nodes,1);

% Automatically set dmesh as a field
params.dmesh = dmesh;

%% Precalculate matrices for least square solvers
params.dmesh = precalc_lsq_matrices(dmesh);

%% Output fields
% This is the list of fields to save in the output files in addition to the
% state variables (hs, zs, hc, Hc, phic, Vm)
params.output_fields = {'phis_x', 'phis_y',...  % Sheet potential gradient
    'dhsdt', 'dhcdt', 'dHcdt', 'dhcdt_area', 'dhsdt_div', 'dhcdt_div',...   % Time derivatives
    'qx_sheet', 'qy_sheet', 'qc',...    % Fluxes
    'dphic_ds', 'phic_node', 'Xi_c',... % Channel potential gradient
    'exchange_frac', 'm_moulin'};       % Moulin input rates
