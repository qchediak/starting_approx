% 2024 Dec 26
% solve the two body system as an ODE

clear; close all

init

% Parameters for G and DG.  Using a cell array structure
% for consistency with other code.
G_params = {d};

% G is the RHS of the ODE y'=G(y), which is the
% two body system, and DG is its Jacobian.
% Here y=(q,p), where q=(q1,q2) and p=(p1,p2).
G = @(y) computeG(y, G_params);
DG = @(y) computeDG(y, G_params);

IRK_params = {h, N, y0, A, b, c, tol};

yns = solve_IRK(G, DG, IRK_params);

plot_qs = true;
plot_ps = true;
plot_energy = true;

plot_title_str = 'Position Coordinates of the Two-Body Problem';
plots_params = {h, N, s, e, alpha0, plot_title_str, method_str, plot_qs, plot_ps, plot_energy};
make_plots(yns, plots_params)

