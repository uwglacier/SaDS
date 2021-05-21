function status = odeID_OutputFcn(t, y, flag, tt, t_wall_ref, its, err)
% odeID_OutputFcn displays a progress bar for integrating an ODE.
%
% status = odeID_OutputFcn(t, y, flag, tt, t_wall_ref, convergence)
%
% The extra argument compared to the matlab syntax lets this function
% display information about the convergence, which is important for an
% explicit solver
%
% See also odeID

if isempty(flag)
    width = length(tt);
    progress = floor(width*(t(1) - tt(1))/(tt(end) - tt(1)));
    pstr = [repmat('.', 1, progress), repmat(' ', 1, width - progress)];
    t2 = toc(t_wall_ref);
    fprintf('[%s] (t_model = %f, t_wall = %f, iter = %d, err = %.3e)\n', pstr, t, t2, its, err)
end
status = 0;
