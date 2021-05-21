function [T,Y] = odeFE(odefun, tspan, dt, v0, opts)
% odeFE provides simple forward Euler timestepping.
%
% odeFE is aa simple Euler timestepping scheme built to use the same
% syntax as the built-in matlab ode integrators (e.g. ode45), but with an
% additional argument to provide the timestep.
%
% [T, Y] = odeRK(odefun, tspan, v0, opts) integrates the ODE odefun for
% the times given by tspan, initial condition v0, and options opts. opts
% must be a structure with fields:
%
%           dt: double
%    outputFcn: function handle
%
% outputFcn must have the signature status = outputFcn(t, y, flag).
% See also progress_function, odeRKplus, odeRK_OutputFcn

%% Handle function inputs
if ~isfield(opts, 'stats')
    opts.stats = ones(size(v0));
end

% Intialize OutputFcn (for consistency with built-in Matlab solvers, see
% https://www.mathworks.com/help/matlab/ref/odeset.html)
if ~isempty(opts.OutputFcn)
    status = opts.OutputFcn(0, v0, 'init', 0);
end

% Expand tspan if it only gives start and end times
if length(tspan)==2
    tspan = tspan(1):opts.dt:tspan(end);
end

%% Time stepping scheme
Y = zeros(length(tspan), length(v0));
T = tspan';
v = v0;
t = tspan(1);
jj=1; % Index for stepping through tspan array
tic;
while t < tspan(end)
    % Determine next time step
    if t+opts.dt >= tspan(jj)
        dt_ii=tspan(jj)-t;
        save_state=true;
    else
        dt_ii=opts.dt;
        save_state=false;
    end

    % This is the advancement of the state
    vprime = odefun(t, v);
    r = max(abs(opts.dt*vprime.*opts.stats));
    v = v + dt_ii*vprime;
    t=t+dt_ii;

    if save_state
        if ~isempty(opts.OutputFcn)
            status = opts.OutputFcn(t, v, [], r);
        end

        Y(jj,:) = v';
        T(jj)=t;
        jj=jj+1;

    end
end

%% Cleanup
if ~isempty(opts.OutputFcn)
    status = opts.OutputFcn([], [], 'done', 0);
end

wall_time = toc;
fprintf('Elapsed time %s seconds\n', wall_time);
end
