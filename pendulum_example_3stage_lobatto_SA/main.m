% 2024 Dec 13
% Solves the pendulum problem as an index 3 DAE using two
% steps with the 3-stage lobatto method and finding the
% error in the starting approximations.
%
% With modifications this could be adapted to code which
% solves a general index 3 problem with n=m=2 and p=1 by a 
% PRK method with s=3.
%
% The purpose of this code is to solve the pendulum problem.  

clear; close all

plots=false;

init 

% Set ratio between step sizes
r=1;

% put everything in a cell array called params
params = {ny, nz, nu, yn, zn, g0, l, A, Ahat, b, c, U10_trivial, U20_trivial, U30_trivial, tol, r};

H_len = 8;
H = logspace(-3.25,-2,H_len);
H = flip(H);

for i=1:H_len
	h = H(i);

	result = find_error_SA(h,params);

	Y0_err = result{1};
	Z0_err = result{2};
	znp0_err = result{3};
	U10_err = result{4};
	U20_err = result{5};
	U30_err = result{6};
	yns = result{7};
	zns = result{8};
	energy_vec = result{9};

	if i==1
		Y0_errs = Y0_err;
		Z0_errs = Z0_err;
		znp0_errs = znp0_err;
		U10_errs = U10_err;
		U20_errs = U20_err;
		U30_errs = U30_err;
	else
		Y0_errs = [Y0_errs; Y0_err];
		Z0_errs = [Z0_errs; Z0_err];
		znp0_errs = [znp0_errs; znp0_err];
		U10_errs = [U10_errs; U10_err];
		U20_errs = [U20_errs; U20_err];
		U30_errs = [U30_errs; U30_err];
	end

	%% Plots 
	if plots && i==H_len
		plots_params = {r,h};
		make_plots(yns,zns,energy_vec,plots_params)
	end
end

Y0_title_str = '$\log_{10} \| Y_{n+1}^{(0)} - Y_{n+1} \|$ vs.~$\log_{10} h$, $r=$';
Z0_title_str = '$\log_{10} \| Z_{n+1}^{(0)} - Z_{n+1} \|$ vs.~$\log_{10} h$, $r=$';
znp0_title_str = '$\log_{10} \| z_{n+1}^{(0)} - z_{n+1} \|$ vs.~$\log_{10} h$, $r=$';
U10_title_str = '$\log_{10} \| U_{n+1,1}^{(0)} - U_{n+1,1} \|$ vs.~$\log_{10} h$, $r=$';
U20_title_str = '$\log_{10} \| U_{n+1,2}^{(0)} - U_{n+1,2} \|$ vs.~$\log_{10} h$, $r=$';
U30_title_str = '$\log_{10} \| U_{n+1,3}^{(0)} - U_{n+1,3} \|$ vs.~$\log_{10} h$, $r=$';

Y0_title_str = append(Y0_title_str, sprintf('%0.1f',r));
Z0_title_str = append(Z0_title_str, sprintf('%0.1f',r));
znp0_title_str = append(znp0_title_str, sprintf('%0.1f',r));
U10_title_str = append(U10_title_str, sprintf('%0.1f',r));
U20_title_str = append(U20_title_str, sprintf('%0.1f',r));
U30_title_str = append(U30_title_str, sprintf('%0.1f',r));

plot_errs_params_Y0 = {Y0_title_str, 'Y0'};
plot_errs_params_Z0 = {Z0_title_str, 'Z0'};
plot_errs_params_znp0 = {znp0_title_str, 'znp0'};
plot_errs_params_U10 = {U10_title_str, 'U10'};
plot_errs_params_U20 = {U20_title_str, 'U20'};
plot_errs_params_U30 = {U30_title_str, 'U30'};

[alpha_Y0, beta_Y0] = plot_initializer_errors(H, Y0_errs, plot_errs_params_Y0);
[alpha_Z0, beta_Z0] = plot_initializer_errors(H, Z0_errs, plot_errs_params_Z0);
[alpha_znp0, beta_znp0] = plot_initializer_errors(H, znp0_errs, plot_errs_params_znp0);
[alpha_U10, beta_U10] = plot_initializer_errors(H, U10_errs, plot_errs_params_U10);
[alpha_U20, beta_U20] = plot_initializer_errors(H, U20_errs, plot_errs_params_U20);
[alpha_U30, beta_U30] = plot_initializer_errors(H, U30_errs, plot_errs_params_U30);

