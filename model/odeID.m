function [T,Y] = odeID(odefun, tspan, v0, opts)
% odeID provides implicit damped timestepping.
%
% odeID is a dampled implicit timestepping scheme built to use the same
% syntax as the build-in matlab ode integrators (e.g. ode45), but with an
% additional argument (a structure that provides information about step
% sizes and tilerances).
%
% [T, Y] = odeID(odefun, tspan, v0, opts) integrates the ODE odefun for
% the times given by tspan, initial condition v0, and options opts. opts
% must be a structure with fields:
%
%           dt: double
%         dtau: double
%      maxIter: Int
%          tol: double
%         visc: double [0 <= visc <= 1]
%    outputFcn: function handle [optional]
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
    status = opts.OutputFcn(0, v0, 'init', 0, 0);
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

% Initialize loop variables
res = zeros(size(v));
v_old = v;
dvdtau = zeros(size(v));
ittot=0;
% Physical time loop
while t < tspan(end)
%     % Determine next time step
    % if t+opts.dt >= tspan(jj)
    %     dt_ii=tspan(jj)-t;
    %     save_state=true;
    % else
    %     dt_ii=opts.dt;
    %     save_state=false;
    % end
    dt_ii = opts.dt;
    save_state = true;
    iter = 0;
    err = 2*opts.tol;
    % Pseudo-transient iteration
        while err>opts.tol && iter<opts.maxIter
            vprime = odefun(t, v);
            res = -(v - v_old)/dt_ii + vprime; % ODE residual

            dvdtau = res + opts.visc*dvdtau;
            v = v + opts.dtau*dvdtau;   % Naive forward Euler timestep
            iter = iter + 1;
            err = norm(res)/length(v);
        end
%     end
    ittot = ittot + iter;
    v_old = v;
    t = t + dt_ii;
    if save_state
        if ~isempty(opts.OutputFcn)
            status = opts.OutputFcn(t, v, [], iter, err);
        end

        Y(jj,:) = v';
        T(jj)=t;
        jj=jj+1;

    end
end

Y(end, :) = v';
%% Cleanup
if ~isempty(opts.OutputFcn)
    status = opts.OutputFcn([], [], 'done', 0, 0);
end

wall_time = toc;
fprintf('Elapsed time %s seconds\n', wall_time);
fprintf('Total iterations: %d\n', ittot)

end
