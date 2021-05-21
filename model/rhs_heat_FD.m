function vprime = rhs_heat_FD(t, v, D, dx)
% rhs_heat_FD. One-dimensional FD heat equation discretization.
%
% vprime = rhs_heat_FD(t, v, D, dx)

q = -D*(v(2:end) - v(1:end-1))/dx;
vprime = -(q(2:end) - q(1:end-1))/dx;
vprime = [0; vprime; 0];