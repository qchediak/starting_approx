% 2025 June 24
% solve the plane pendulum as a DAE and plot starting approximation
% errors

clear; close all

init

r = 1;

ynm = y0;
znm = z0;
unm = u0;

x0_trivial = [kron(ones(s,1),ynm); kron(ones(s,1),znm); kron(ones(s,1),unm)];

% parameters for find_error_SA
params = {ny, nz, nu, g0, l, b, c, A, Ahat, ynm, znm, unm, x0_trivial, tol, r, method_str, f, k, g};

stepsizes_len = 8;
stepsizes = logspace(-2.25,-1.5,stepsizes_len);
stepsizes = flip(stepsizes);

for i=1:stepsizes_len
	h = stepsizes(i);

	result = find_error_SA(h, params);

	% unpack result
	Ynp0_err = result{1};
	Znp0_err = result{2};
	Unp0_err = result{3};
	yns = result{4};
	zns = result{5};
	energy_vec = result{6};

	if i==1
		Y0_errs = Ynp0_err;
		Z0_errs = Znp0_err;
		U0_errs = Unp0_err;
	else
		Y0_errs = [Y0_errs; Ynp0_err];
		Z0_errs = [Z0_errs; Znp0_err];
		U0_errs = [U0_errs; Unp0_err];
	end
end

% method_str_edit is method_str in lowercase with no spaces
method_str_edit = lower(method_str);
method_str_edit = strrep(method_str_edit, ' ', '');

line1_Y = '$\log_{10} \| Y_{n+1}^{(0)} - Y_{n+1} \|$ vs.~$\log_{10} h$';
line1_Z = '$\log_{10} \| Z_{n+1}^{(0)} - Z_{n+1} \|$ vs.~$\log_{10} h$';
line1_U = '$\log_{10} \| U_{n+1}^{(0)} - U_{n+1} \|$ vs.~$\log_{10} h$';

line2 =  sprintf('%s, $s=%i$, $r=%.1f$', method_str, s, r);

title_str_Y = sprintf('%s \n %s', line1_Y, line2);
title_str_Z = sprintf('%s \n %s', line1_Z, line2);
title_str_U = sprintf('%s \n %s', line1_U, line2);

print_str_Y = 'Y0';
print_str_Z = 'Z0';
print_str_U = 'U0';

plot_errs_params_Y = {title_str_Y, print_str_Y};
plot_errs_params_Z = {title_str_Z, print_str_Z};
plot_errs_params_U = {title_str_U, print_str_U};

[alpha_Y0, beta_Y0] = plot_initializer_errors(stepsizes, Y0_errs, plot_errs_params_Y);
[alpha_Z0, beta_Z0] = plot_initializer_errors(stepsizes, Z0_errs, plot_errs_params_Z);
[alpha_U0, beta_U0] = plot_initializer_errors(stepsizes, U0_errs, plot_errs_params_U);


