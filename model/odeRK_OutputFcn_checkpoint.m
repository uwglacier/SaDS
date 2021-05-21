function status = odeRK_OutputFcn_checkpoint(t, y, flag, tt, checkpoints, t_wall_ref, convergence)
% odeRK_OutputFcn_checkpoint displays a progress bar for integrating an ODE
% while saving outputs at checkpoint times
%
% status = progress_function(t, y, flag, tt, t_wall_ref, convergence)
%
% The extra argument compared to the matlab syntax lets this function
% display information about the convergence, which is important for an
% explicit solver
%
% See also odeRK

if isempty(flag)
    width = length(tt);
    progress = floor(width*(t(1) - tt(1))/(tt(end) - tt(1)));
    pstr = [repmat('.', 1, progress), repmat(' ', 1, width - progress)];
    t2 = toc(t_wall_ref);
    fprintf('[%s] (t_model = %f, t_wall = %f, AbsConv = %f)\n', pstr, t, t2, convergence)
    
    if ismember(t, checkpoints)
        jj = find(checkpoints==t);
        save(sprintf('checkpoint_%03d.mat', jj))
    end
        
end
status = 0;