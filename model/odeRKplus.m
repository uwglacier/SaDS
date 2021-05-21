function [T,Y] = odeRKplus(odefun, tspan, v0, opts)
% odeRKplus provides naive adaptive Runge-Kutta timestepping.
%
% [T, Y] = odeRKplus(odefun, tspan, v0, opts) integrates the ODE odefun for
% the times given by tspan, initial condition v0, and options opts. opts
% must be a structure with fields:
%
%           dt: double
%       dt_min: double
%       dt_max: double
%       AbsTol: double
%    OutputFcn: function handle
%
% odeRKplus calculates the timestep required to reach a tolerance of AbsTol
% using the previous timestep's derivative, and enforcing dt_min <= dt <=
% dt_max. OutputFcn has the same usage as in built-in odeXX functions.
%
% OutputFcn must have the signature status = outputFcn(t, y, flag).
% See also progress_function, odeRK

%% Handle function inputs
% Intialize OutputFcn (for consistency with built-in Matlab solvers, see
% https://www.mathworks.com/help/matlab/ref/odeset.html)
if ~isempty(opts.OutputFcn)
    status = opts.OutputFcn(0, v0, 'init');
end

% Unpack opts
dt_min = opts.dt_min;
dt_max = opts.dt_max;
r = opts.AbsTol;

% Define how we calculate the timestep
get_dt = @(vprime) min( max(dt_min, min(abs(r./vprime))), dt_max);

%% Time stepping scheme
t=tspan(1);
dt = opts.dt_min;
v = v0;
Y=zeros(length(tspan),length(v0));
T=tspan';
jj=1; % Index for stepping through tspan array
n_steps = 0;    % Keep track of the total steps taken
tic;
while t < tspan(end)
    % Determine next time step
    
    % Case 1 - we save outputs at specific times
    if length(tspan)>2
        if t+dt >= tspan(jj)
            dt_ii=tspan(jj)-t;
            save_state = true;
        else
            dt_ii=dt;
            save_state = false;
        end
    else
        dt_ii = dt;
        save_state = true;
    end
    % This is the advancement of the state
    k1 = odefun(t, v);
    k2 = odefun(t + dt_ii/2, v + dt_ii*k1/2);
    k3 = odefun(t + dt_ii/2, v + dt_ii*k2/2);
    k4 = odefun(t + dt_ii, v + dt_ii*k3);
    
    vprime = (k1 + 2*k2 + 2*k3 + k4)/6;
    
    v = v + dt_ii*vprime;
    t=t+dt_ii;
    
    % Recalculate dt and update step counter
    dt = get_dt(vprime.*opts.stats);
    n_steps = n_steps+1;
    
    if save_state
        if ~isempty(opts.OutputFcn)
            status = opts.OutputFcn(t, v, []);
        end
        Y(jj,:) = v';
        T(jj)=t;
        jj=jj+1;
    end
end

%% Clean up
% Cleanup tasks for OutputFcn
if ~isempty(opts.OutputFcn)
    status = opts.OutputFcn([], [], 'done');
end

toc
end
