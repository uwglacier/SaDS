function status = progress_function(t, y, flag, tt, t_wall_ref)
% progress_function displays a progress bar for integrating an ODE.
%
% status = progress_function(t, y, flag, tt, t_wall_ref)
%
% Use as follows with built-in matlab solvers:
%
% t0 = tic;
% Fcn = @(t, y, flag) progress_function(t, y, flag, tspan, t0);
% odeset('OutputFcn', Fcn);

if isempty(flag)
    width = length(tt);
    progress = floor(width*t/tt(end));
    pstr = [repmat('.', 1, progress), repmat(' ', 1, width - progress)];
    t2 = toc(t_wall_ref);
    fprintf('[%s] (t_model = %f, t_wall = %f)\n', pstr, t, t2)
end
status = 0;