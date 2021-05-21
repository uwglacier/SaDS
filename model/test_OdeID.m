% Script to test damped iterative ODE solver
%
% Based off of EGU presentation. This is a finite difference
% one-dimensional heat equation code

lx = 10;
D = 1;
ttot = 0.6;
dt = 0.1;

nx = 256;
tol = 1e-6;
maxIter = 1e5;
damp = 1 - 41/nx;

dx = lx/nx;
dtau = (1.0/(dx^2/D/2.1) + 1.0/dt)^-1;
xc = dx/2:dx:(lx-dx/2);

dt_cfl = dx^2/D/2.1;

H0 = exp(-(xc - lx/2).^2)';

tt =0:dt:ttot;


odefun = @(t, Y) rhs_heat_FD(t, Y, D, dx);

solver_opts.dt = dt_cfl;
solver_opts.OutputFcn = [];
[T, Y] = odeRK(odefun, tt, H0, solver_opts);

figure
hold on
plot(xc, Y')
v_analytic = 1/sqrt(4*(ttot+1/4)) * exp(-(xc-lx/2).^2 /(4*(ttot+1/4)));
err = norm(Y(end, :) - v_analytic)

opts_damp.dt = 0.1;
opts_damp.dtau = dtau;
opts_damp.maxIter = maxIter;
opts_damp.visc = damp;
opts_damp.tol = tol;

tref = tic;
opts_damp.OutputFcn = @(t, y, flag, ittot, err) odeID_OutputFcn(t, y, flag, tt, tref, ittot, err);
[T, Y_ID] = odeID(odefun, tt, H0, opts_damp);

figure
hold on
plot(xc, Y_ID')
% plot(xc, H0, 'r')
plot(xc, v_analytic, 'k--')

figure
hold on
plot(xc, v_analytic - Y_ID(end, :))
err = norm(Y(end, :) - v_analytic)
