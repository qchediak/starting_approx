% 2024 Dec 24
% solve the plane pendulum as an ODE and plot starting approximation
% errors

clear; close all

init

r = 1/2;

% parameters for G and DG
G_params = {g0, ny, nz, nu}; 

G = @(psi) computeG(psi, G_params);
DG = @(psi) computeDG(psi, G_params);

psinm = psi0; % initial value

% parameters for find_error_SA
params = {ny, nz, nu, g0, l, A, b, c, psinm, Psi0_trivial, tol, r, method_str};

H_len = 10;
H = logspace(-4,-2,H_len);
H = flip(H);

for i=1:H_len
	h = H(i);

	result = find_error_SA(h, params);

	Psinp0_err = result{1};
	psins = result{2};
	energy_vec = result{3};

	if i==1
		Psi0_errs = Psinp0_err;
	else
		Psi0_errs = [Psi0_errs; Psinp0_err];
	end
end

% method_str_edit is method_str in lowercase with no spaces
method_str_edit = lower(method_str);
method_str_edit = strrep(method_str_edit, ' ', '');

if strcmp(method_str_edit, lower('lobattoIIIA')) && s==2
	line1 = '$\log_{10} \| \Psi_{n+1,2}^{(0)} - \Psi_{n+1,2} \|$ vs.~$\log_{10} h$';
else
	line1 = '$\log_{10} \| \Psi_{n+1}^{(0)} - \Psi_{n+1} \|$ vs.~$\log_{10} h$';
end

line2 =  sprintf('%s, $s=%i$, $r=%.1f$', method_str, s, r);
title_str = sprintf('%s \n %s', line1, line2);
plot_errs_params = {title_str, 'Psi0'};

[alpha_Y0, beta_Y0] = plot_initializer_errors(H, Psi0_errs, plot_errs_params);

