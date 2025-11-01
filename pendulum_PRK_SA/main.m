% 2025 June 24
% solve the plane pendulum as a DAE and plot starting approximation
% errors

clear; close all

init

r = 1/2;

ynm = y0;
znm = z0;
unm = u0;

x0_trivial = [kron(ones(s,1),ynm); kron(ones(s,1),znm); kron(ones(s,1),unm)];

% parameters for find_error_SA
params = {ny, nz, nu, g0, l, b, c, A, Ahat, ynm, znm, unm, x0_trivial, tol, r, method_str, f, k, g};

stepsizes_len = 8;
stepsizes = logspace(-3.3,-3,stepsizes_len);
stepsizes = flip(stepsizes);

for i=1:stepsizes_len
	h = stepsizes(i);

	result = find_error_SA(h, params);

	% unpack result
	Ynp0_err_vec = result{1};
	Znp0_err_vec = result{2};
	Unp0_err_vec = result{3};
	yns = result{4};
	zns = result{5};
	energy_vec = result{6};

	% Keep errors in either component or composite form depending on
	% error plotting option set in init file.
	% If the option is set to 'composite', take the norm of the entire
	% error vector.  Otherwise, take the norm of each individual
	% component.

	% Y component
	if strcmpi(error_plotting_Y, 'composite')
		Ynp0_err = norm(Ynp0_err_vec);
		if i==1
			Y0_errs = Ynp0_err;
		else
			Y0_errs = [Y0_errs; Ynp0_err];
		end
	elseif strcmpi(error_plotting_Y, 'component') && s==3
		Ynp10_err = norm(Ynp0_err_vec(1:ny));			% if ny=2, this is 1:2
		Ynp20_err = norm(Ynp0_err_vec(ny+1:2*ny));		% if ny=2, this is 3:4
		Ynp30_err = norm(Ynp0_err_vec(2*ny+1:3*ny));	% if ny=2, this is 5:6
		if i==1
			Y10_errs = Ynp10_err;
			Y20_errs = Ynp20_err;
			Y30_errs = Ynp30_err;
		else
			Y10_errs = [Y10_errs; Ynp10_err];
			Y20_errs = [Y20_errs; Ynp20_err];
			Y30_errs = [Y30_errs; Ynp30_err];
		end
	else
		error('Error (main): Invalid combination of error_plotting_Y and s')
	end

	% Z component
	if strcmpi(error_plotting_Z, 'composite')
		Znp0_err = norm(Znp0_err_vec);
		if i==1
			Z0_errs = Znp0_err;
		else
			Z0_errs = [Z0_errs; Znp0_err];
		end
	elseif strcmpi(error_plotting_Z, 'component') && s==3
		Znp10_err = norm(Znp0_err_vec(1:nz));			% if nz=2, this is 1:2
		Znp20_err = norm(Znp0_err_vec(nz+1:2*nz));		% if nz=2, this is 3:4
		Znp30_err = norm(Znp0_err_vec(2*nz+1:3*nz));	% if nz=2, this is 5:6
		if i==1
			Z10_errs = Znp10_err;
			Z20_errs = Znp20_err;
			Z30_errs = Znp30_err;
		else
			Z10_errs = [Z10_errs; Znp10_err];
			Z20_errs = [Z20_errs; Znp20_err];
			Z30_errs = [Z30_errs; Znp30_err];
		end
	else
		error('Error (main): Invalid combination of error_plotting_Z and s')
	end

	% U component
	if strcmpi(error_plotting_U, 'composite')
		Unp0_err = norm(Unp0_err_vec);
		if i==1
			U0_errs = Unp0_err;
		else
			U0_errs = [U0_errs; Unp0_err];
		end
	elseif strcmpi(error_plotting_U, 'component') && s==3
		Unp10_err = norm(Unp0_err_vec(1:nu));			% if nu=1, this is 1:1
		Unp20_err = norm(Unp0_err_vec(nu+1:2*nu));		% if nu=1, this is 2:2
		Unp30_err = norm(Unp0_err_vec(2*nu+1:3*nu));	% if nu=1, this is 3:3
		if i==1
			U10_errs = Unp10_err;
			U20_errs = Unp20_err;
			U30_errs = Unp30_err;
		else
			U10_errs = [U10_errs; Unp10_err];
			U20_errs = [U20_errs; Unp20_err];
			U30_errs = [U30_errs; Unp30_err];
		end
	else
		error('Error (main): Invalid combination of error_plotting_U and s')
	end
end

% method_str_edit is method_str in lowercase with no spaces
method_str_edit = lower(method_str);
method_str_edit = strrep(method_str_edit, ' ', '');

line2 =  sprintf('%s, $s=%i$, $r=%.1f$', method_str, s, r);

% Y component
if strcmpi(error_plotting_Y, 'composite')
	line1_Y = '$\log_{10} \| Y_{n+1}^{(0)} - Y_{n+1} \|$ vs.~$\log_{10} h$';
	title_str_Y = sprintf('%s \n %s', line1_Y, line2);
	print_str_Y = 'Y0';
	plot_errs_params_Y = {title_str_Y, print_str_Y};
	[alpha_Y0, beta_Y0] = plot_initializer_errors(stepsizes, Y0_errs, plot_errs_params_Y);
elseif strcmpi(error_plotting_Y, 'component') && s==3
	line1_Y1 = '$\log_{10} \| Y_{n+1,1}^{(0)} - Y_{n+1,1} \|$ vs.~$\log_{10} h$';
	title_str_Y1 = sprintf('%s \n %s', line1_Y1, line2);
	print_str_Y1 = 'Y10';
	plot_errs_params_Y1 = {title_str_Y1, print_str_Y1};
	[alpha_Y10, beta_Y10] = plot_initializer_errors(stepsizes, Y10_errs, plot_errs_params_Y1);

	line1_Y2 = '$\log_{10} \| Y_{n+1,2}^{(0)} - Y_{n+1,2} \|$ vs.~$\log_{10} h$';
	title_str_Y2 = sprintf('%s \n %s', line1_Y2, line2);
	print_str_Y2 = 'Y20';
	plot_errs_params_Y2 = {title_str_Y2, print_str_Y2};
	[alpha_Y20, beta_Y20] = plot_initializer_errors(stepsizes, Y20_errs, plot_errs_params_Y2);

	line1_Y3 = '$\log_{10} \| Y_{n+1,3}^{(0)} - Y_{n+1,3} \|$ vs.~$\log_{10} h$';
	title_str_Y3 = sprintf('%s \n %s', line1_Y3, line2);
	print_str_Y3 = 'Y30';
	plot_errs_params_Y3 = {title_str_Y3, print_str_Y3};
	[alpha_Y30, beta_Y30] = plot_initializer_errors(stepsizes, Y30_errs, plot_errs_params_Y3);
else
	error('Error (main): Invalid combination of error_plotting_Y and s')
end

% Z component
if strcmpi(error_plotting_Z, 'composite')
	line1_Z = '$\log_{10} \| Z_{n+1}^{(0)} - Z_{n+1} \|$ vs.~$\log_{10} h$';
	title_str_Z = sprintf('%s \n %s', line1_Z, line2);
	print_str_Z = 'Z0';
	plot_errs_params_Z = {title_str_Z, print_str_Z};
	[alpha_Z0, beta_Z0] = plot_initializer_errors(stepsizes, Z0_errs, plot_errs_params_Z);
elseif strcmpi(error_plotting_Z, 'component') && s==3
	line1_Z1 = '$\log_{10} \| Z_{n+1,1}^{(0)} - Z_{n+1,1} \|$ vs.~$\log_{10} h$';
	title_str_Z1 = sprintf('%s \n %s', line1_Z1, line2);
	print_str_Z1 = 'Z10';
	plot_errs_params_Z1 = {title_str_Z1, print_str_Z1};
	[alpha_Z10, beta_Z10] = plot_initializer_errors(stepsizes, Z10_errs, plot_errs_params_Z1);

	line1_Z2 = '$\log_{10} \| Z_{n+1,2}^{(0)} - Z_{n+1,2} \|$ vs.~$\log_{10} h$';
	title_str_Z2 = sprintf('%s \n %s', line1_Z2, line2);
	print_str_Z2 = 'Z20';
	plot_errs_params_Z2 = {title_str_Z2, print_str_Z2};
	[alpha_Z20, beta_Z20] = plot_initializer_errors(stepsizes, Z20_errs, plot_errs_params_Z2);

	line1_Z3 = '$\log_{10} \| Z_{n+1,3}^{(0)} - Z_{n+1,3} \|$ vs.~$\log_{10} h$';
	title_str_Z3 = sprintf('%s \n %s', line1_Z3, line2);
	print_str_Z3 = 'Z30';
	plot_errs_params_Z3 = {title_str_Z3, print_str_Z3};
	[alpha_Z30, beta_Z30] = plot_initializer_errors(stepsizes, Z30_errs, plot_errs_params_Z3);
else
	error('Error (main): Invalid combination of error_plotting_Z and s')
end

% U component
if strcmpi(error_plotting_U, 'composite')
	line1_U = '$\log_{10} \| U_{n+1}^{(0)} - U_{n+1} \|$ vs.~$\log_{10} h$';
	title_str_U = sprintf('%s \n %s', line1_U, line2);
	print_str_U = 'U0';
	plot_errs_params_U = {title_str_U, print_str_U};
	[alpha_U0, beta_U0] = plot_initializer_errors(stepsizes, U0_errs, plot_errs_params_U);
elseif strcmpi(error_plotting_U, 'component') && s==3
	line1_U1 = '$\log_{10} \| U_{n+1,1}^{(0)} - U_{n+1,1} \|$ vs.~$\log_{10} h$';
	title_str_U1 = sprintf('%s \n %s', line1_U1, line2);
	print_str_U1 = 'U10';
	plot_errs_params_U1 = {title_str_U1, print_str_U1};
	[alpha_U10, beta_U10] = plot_initializer_errors(stepsizes, U10_errs, plot_errs_params_U1);

	line1_U2 = '$\log_{10} \| U_{n+1,2}^{(0)} - U_{n+1,2} \|$ vs.~$\log_{10} h$';
	title_str_U2 = sprintf('%s \n %s', line1_U2, line2);
	print_str_U2 = 'U20';
	plot_errs_params_U2 = {title_str_U2, print_str_U2};
	[alpha_U20, beta_U20] = plot_initializer_errors(stepsizes, U20_errs, plot_errs_params_U2);

	line1_U3 = '$\log_{10} \| U_{n+1,3}^{(0)} - U_{n+1,3} \|$ vs.~$\log_{10} h$';
	title_str_U3 = sprintf('%s \n %s', line1_U3, line2);
	print_str_U3 = 'U30';
	plot_errs_params_U3 = {title_str_U3, print_str_U3};
	[alpha_U30, beta_U30] = plot_initializer_errors(stepsizes, U30_errs, plot_errs_params_U3);
else
	error('Error (main): Invalid combination of error_plotting_U and s')
end




