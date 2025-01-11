% 2024 Dec 28
% solve the two body problem as an ODE and plot starting approximation
% errors

clear; close all

init

r = 1/2;

% Parameters for G and DG.  Using a cell array structure
% for consistency with other code.
G_params = {d};

% G is the RHS of the ODE y'=G(y), which is the
% two body system, and DG is its Jacobian.
% Here y=(q,p), where q=(q1,q2) and p=(p1,p2).
G = @(y) computeG(y, G_params);
DG = @(y) computeDG(y, G_params);

ynm = y0; % initial value

% parameters for find_error_SA
params = {d, A, b, c, ynm, Y0_trivial, tol, r, method_str}; 

H_len = 10;
H = logspace(-4,-2,H_len);
H = flip(H);

for i=1:H_len
	h = H(i);

	result = find_error_SA(h, params);

	Ynp0_err = result{1};
	yns = result{2};
	energy_vec = result{3};

	if i==1
		Y0_errs = Ynp0_err;
	else
		Y0_errs = [Y0_errs; Ynp0_err];
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
plot_errs_params = {title_str, 'Y0'};

[alpha_Y0, beta_Y0] = plot_initializer_errors(H, Y0_errs, plot_errs_params);


