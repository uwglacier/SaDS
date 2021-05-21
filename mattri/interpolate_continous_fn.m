function value_nodes2 = interpolate_continous_fn(mesh1, mesh2, value_nodes1)
% value_nodes2 = interpolate_continous_fn(mesh1, mesh2, value_nodes1)
%
% Takes a mesh1 and associated function values value_nodes1 as well as
% another mesh2. Returns the interpolated function values on the new mesh

% This m-file is part of the mattri meshing package.
% Copyright (c) 2014, Mauro A Werder
% All rights reserved.
% Licensed under a BSD 2-Clause License, see LICENCE file

interp_fn = TriScatteredInterp(mesh1.nodes, value_nodes1);

value_nodes2 = interp_fn(mesh2.nodes);