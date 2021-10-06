function model_outputs = run_model(params, Y0, varargin)
% run_model wraps the core model code to provide a more user-friendly way
% to run the model.
% Inputs:
%   * params: parameter structure (e.g. from get_default_params)
%   * Y0: structure providing initial conditions on hs, zs, hc, Hc, phic
%   * outputs [optional]: path to save model outputs. If not passed,
%                           outputs are not saved to file
% Returns:
%   * model_outputs: structure with fields
%       params: copy of input params
%       outputs: structure containing model outputs, with fields
%                hs: [n_elements × n_outputs double]
%                zs: [n_elements × n_outputs double]
%                hc: [n_edges × n_outputs double]
%                Hc: [n_edges × n_outputs double]
%                zc: [n_edges × n_outputs double]
%                tt: [1 × n_outputs double]
%                Vm: [n_moulins × n_outputs double]
%            phis_x: [n_elements × n_outputs double]
%            phis_y: [n_elements × n_outputs double]
%        phis_bndry: [n_edges × n_outputs double]
%             dhsdt: [n_elements × n_outputs double]
%             dhcdt: [n_edges × n_outputs double]
%             dHcdt: [n_edges × n_outputs double]
%         dhsdt_div: [n_elements × n_outputs double]
%         dhcdt_div: [n_edges × n_outputs double]
%           qx_edge: [n_edges × n_outputs double]
%           qy_edge: [n_edges × n_outputs double]
%          qx_sheet: [n_elements × n_outputs double]
%          qy_sheet: [n_elements × n_outputs double]
%                qc: [n_edges × n_outputs double]
%              phic: [n_edges × n_outputs double]
%          dphic_ds: [n_edges × n_outputs double]
%              Xi_c: [n_edges × n_outputs double]
%     exchange_frac: [n_edges × n_outputs double]
%          m_moulin: [n_moulins × n_outputs double]

if ~isempty(varargin)
    filename = varargin{1};
else
    filename = 'none';
end

% Make sure that the inputs are function handles (this ensures backwards
% compatability)
if ~isa(params.ms, 'function_handle')
    params.ms = @(t) params.ms;
end

if ~isa(params.msc, 'function_handle')
    params.msc = @(t) params.msc;
end

%% Validate the params
status = validate_params(params);

% Print parameter structure
params

if params.stop_on_warning && status==0
    error('Parameter validation failed')
end

%% Construct ode function -------------------------------------------------
odefun = @(t,y) rhs_sheet_and_channel(t, y, params.dmesh, params, 'solver');

% Initial condition state vector
v0 = [Y0.hs; Y0.zs; Y0.hc; Y0.Hc; Y0.phic; Y0.Vm];

% Define the output function (print progress to terminal)
t_wall_ref = tic;
OutFcn = @(t, y, flag) progress_function(t, y, flag, params.tt, t_wall_ref);
stats = [ones(size(Y0.hs)); 0*Y0.zs;...
        ones(size(Y0.hc)); 0*Y0.Hc; 0*Y0.phic; 0*Y0.Vm];

if strcmp(params.solver,'odeRK')
    disp('Using RK time stepping')

    % Add the appropriate output function to display solver progress
    params.solver_opts.OutputFcn = @(t, y, flag, r) odeRK_OutputFcn(t, y, flag, params.tt, t_wall_ref, r);
    params.solver_opts.stats = stats;

    [~, Y] = odeRK(odefun, params.tt, v0, params.solver_opts);

elseif strcmp(params.solver,'odeFE')
        disp('Using FE time stepping')

        % Add the appropriate output function to display solver progress
        params.solver_opts.OutputFcn = @(t, y, flag, r) odeRK_OutputFcn(t, y, flag, params.tt, t_wall_ref, r);
        params.solver_opts.stats = stats;

        [~, Y] = odeFE(odefun, params.tt, v0, params.solver_opts);

elseif strcmp(params.solver, 'odeID')
    disp('Using iterative damped time stepping')

    % Add the appropriate output function to display solver progress
    params.solver_opts.OutputFcn = @(t, y, flag, iter, err) odeID_OutputFcn(t, y, flag, params.tt, t_wall_ref, iter, err);

    [~, Y] = odeID(odefun, params.tt, v0, params.solver_opts);

else
    % Attempt to get the ode solver function handle
    odesolver = str2func(params.solver);
    fprintf('Using %s time stepping\n', char(odesolver));

    % Add the appropriate output function to display solver progress
    params.solver_opts.OutputFcn = OutFcn;

    [~, Y] = odesolver(odefun, params.tt, v0, params.solver_opts);

end

% Now save the outputs
model_outputs = save_model_outputs(filename, params, Y);
