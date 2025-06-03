% 2024 Dec 21
% solve the plane pendulum as an ODE

clear; close all

init

h = 0.1;
N = 50;

% parameters for G and DG
G_params = {g0, ny, nz, nu}; 

G = @(psi) computeG(psi, G_params);
DG = @(psi) computeDG(psi, G_params);

%f = @(y) 2*y; % test function
%Df = @(y) 2*eye(n);

psins = solve_IRK(G,DG,h,N,psi0,A,b,c,tol);

plot_ys = true;
plot_zs = false;
plot_energy = true;

plot_title_str = 'Position Coordinates $(y_1,y_2)$ of the Pendulum';
plots_params = {h, N, s, g0, l, plot_title_str, method_str, plot_ys, plot_zs, plot_energy};
make_plots(psins, plots_params)

