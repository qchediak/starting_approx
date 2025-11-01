% 2025 Oct 31
% solve the plane pendulum as a DAE using projected Lobatto IIIC and
% plot starting approximation errors

clear; close all

init

r = 1;

ynm = y0;
znm = z0;
unm = u0;

x0_trivial = [kron(ones(s,1),ynm); kron(ones(s,1),znm); kron(ones(s,1),unm)];

% parameters for find_error_SA
params = {ny, nz, nu, g0, l, b, c, A, Ahat, ynm, znm, unm, x0_trivial, tol, r, method_str, f, k, g, gyf, iterative_scheme1, iterative_scheme2};

stepsizes_len = 8;
stepsizes = logspace(-3.3,-3,stepsizes_len);
%stepsizes = logspace(-3.5,-3.2,stepsizes_len);
stepsizes = flip(stepsizes);

% -----------------------------------------------------
% Initialize the error vectors as vectors of NaN values
% -----------------------------------------------------

if error_plotting_option==1
	if strcmpi(error_plotting_Y, 'composite')
		Y0_errs = NaN(stepsizes_len,1);
	elseif strcmpi(error_plotting_Y, 'component') && s==3
		Y10_errs = NaN(stepsizes_len,1);
		Y20_errs = NaN(stepsizes_len,1);
		Y30_errs = NaN(stepsizes_len,1);
	else
		error('Error (main): Invalid value of error_plotting_Y')
	end

	if strcmpi(error_plotting_Z, 'composite')
		Z0_errs = NaN(stepsizes_len,1);
	elseif strcmpi(error_plotting_Z, 'component') && s==3
		Z10_errs = NaN(stepsizes_len,1);
		Z20_errs = NaN(stepsizes_len,1);
		Z30_errs = NaN(stepsizes_len,1);
	else
		error('Error (main): Invalid value of error_plotting_Z')
	end

	if strcmpi(error_plotting_U, 'composite')
		U0_errs = NaN(stepsizes_len,1);
	elseif strcmpi(error_plotting_U, 'component') && s==3
		U10_errs = NaN(stepsizes_len,1);
		U20_errs = NaN(stepsizes_len,1);
		U30_errs = NaN(stepsizes_len,1);
	else
		error('Error (main): Invalid value of error_plotting_U')
	end

	znp0_errs = NaN(stepsizes_len,1);
	Unpsp0_errs = NaN(stepsizes_len,1);
elseif error_plotting_option==2 && s==3
	Y10_errs = NaN(stepsizes_len,1);
	Y20_errs = NaN(stepsizes_len,1);
	Y30_errs = NaN(stepsizes_len,1);

	Z0_errs = NaN(stepsizes_len,1);
	U0_errs = NaN(stepsizes_len,1);
else
	error('Error (main): Invalid combination of error_plotting_option and s.')
end

% ----------
% Main loop
% ----------

for i=1:stepsizes_len
	h = stepsizes(i);

	result = find_error_SA(h, params);

	% unpack result
	Ynp0_err_vec = result{1};
	Znp0_err_vec = result{2};
	Unp0_err_vec = result{3};
	znp0_err_vec = result{4};
	Unpsp0_err_vec = result{5};
	yns = result{6};
	zns = result{7};
	energy_vec = result{8};

	% update error vector

	if error_plotting_option==1
		% Find vector of errors for znp0
		znp0_err = norm(znp0_err_vec);
		znp0_errs(i) = znp0_err;

		% Find vector of errors for Unpsp0
		Unpsp0_err = norm(Unpsp0_err_vec);
		Unpsp0_errs(i) = Unpsp0_err;

		% Find errors for Y0, Z0, and U0.
		%
		% Keep errors in either component or composite form
		% depending on error plotting option set in init file.
		% If the option is set to 'composite', take the norm of the entire
		% error vector.  Otherwise, take the norm of each individual
		% component.

		% Y component
		if strcmpi(error_plotting_Y, 'composite')
			Ynp0_err = norm(Ynp0_err_vec);
			Y0_errs(i) = Ynp0_err;
		elseif strcmpi(error_plotting_Y, 'component') && s==3
			Ynp10_err = norm(Ynp0_err_vec(1:ny));			% if ny=2, this is 1:2
			Ynp20_err = norm(Ynp0_err_vec(ny+1:2*ny));		% if ny=2, this is 3:4
			Ynp30_err = norm(Ynp0_err_vec(2*ny+1:3*ny));	% if ny=2, this is 5:6

			Y10_errs(i) = Ynp10_err;
			Y20_errs(i) = Ynp20_err;
			Y30_errs(i) = Ynp30_err;
		else
			error('Error (main): Invalid combination of error_plotting_Y and s')
		end

		% Z component
		if strcmpi(error_plotting_Z, 'composite')
			Znp0_err = norm(Znp0_err_vec);
			Z0_errs(i) = Znp0_err;
		elseif strcmpi(error_plotting_Z, 'component') && s==3
			Znp10_err = norm(Znp0_err_vec(1:nz));			% if nz=2, this is 1:2
			Znp20_err = norm(Znp0_err_vec(nz+1:2*nz));		% if nz=2, this is 3:4
			Znp30_err = norm(Znp0_err_vec(2*nz+1:3*nz));	% if nz=2, this is 5:6

			Z10_errs(i) = Znp10_err;
			Z20_errs(i) = Znp20_err;
			Z30_errs(i) = Znp30_err;
		else
			error('Error (main): Invalid combination of error_plotting_Z and s')
		end

		% U component
		if strcmpi(error_plotting_U, 'composite')
			Unp0_err = norm(Unp0_err_vec);
			U0_errs(i) = Unp0_err;
		elseif strcmpi(error_plotting_U, 'component') && s==3
			Unp10_err = norm(Unp0_err_vec(1:nu));			% if nu=1, this is 1:1
			Unp20_err = norm(Unp0_err_vec(nu+1:2*nu));		% if nu=1, this is 2:2
			Unp30_err = norm(Unp0_err_vec(2*nu+1:3*nu));	% if nu=1, this is 3:3

			U10_errs(i) = Unp10_err;
			U20_errs(i) = Unp20_err;
			U30_errs(i) = Unp30_err;
		else
			error('Error (main): Invalid combination of error_plotting_U and s')
		end
	elseif error_plotting_option==2
		% Y component: Keep errors for each component separate
		Ynp10_err = norm(Ynp0_err_vec(1:ny));			% if ny=2, this is 1:2
		Ynp20_err = norm(Ynp0_err_vec(ny+1:2*ny));		% if ny=2, this is 3:4
		Ynp30_err = norm(Ynp0_err_vec(2*ny+1:3*ny));	% if ny=2, this is 5:6

		Y10_errs(i) = Ynp10_err;
		Y20_errs(i) = Ynp20_err;
		Y30_errs(i) = Ynp30_err;

		% Z component: Combine Znp0 with znp0 
		Znp0_err_vec = [Znp0_err_vec; znp0_err_vec];
		Znp0_err = norm(Znp0_err_vec);
		Z0_errs(i) = Znp0_err;

		% U component: Combine Unp0 with Unpsp0 with U0
		Unp0_err_vec = [Unp0_err_vec; Unpsp0_err_vec];
		Unp0_err = norm(Unp0_err_vec);
		U0_errs(i) = Unp0_err;
	else
		error('Error (main): Invalid value of error_plotting_option)')
	end
end

% ----------------
% Plot the errors
% ----------------

% method_str_edit is method_str in lowercase with no spaces
method_str_edit = lower(method_str);
method_str_edit = strrep(method_str_edit, ' ', '');

line2 =  sprintf('%s, $s=%i$, $r=%.1f$', method_str, s, r);

if error_plotting_option==1
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

	% znp

	line1_znp = '$\log_{10} \| z_{n+1}^{(0)} - z_{n+1} \|$ vs.~$\log_{10} h$';
	title_str_znp = sprintf('%s \n %s', line1_znp, line2);
	print_str_znp = 'znp0';
	plot_errs_params_znp = {title_str_znp, print_str_znp};
	[alpha_znp0, beta_znp0] = plot_initializer_errors(stepsizes, znp0_errs, plot_errs_params_znp);


	% Unpsp

	line1_Unpsp = '$\log_{10} \| U_{n+1,s+1}^{(0)} - U_{n+1,s+1} \|$ vs.~$\log_{10} h$';
	title_str_Unpsp = sprintf('%s \n %s', line1_Unpsp, line2);
	print_str_Unpsp = 'Unpsp0';
	plot_errs_params_Unpsp = {title_str_Unpsp, print_str_Unpsp};
	[alpha_Unpsp0, beta_Unpsp0] = plot_initializer_errors(stepsizes, Unpsp0_errs, plot_errs_params_Unpsp);

elseif error_plotting_option==2 && s==3
	% Y component
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

	% Z component
	line1_Z = '$\log_{10} \| Z_{n+1}^{(0)} - Z_{n+1} \|$ vs.~$\log_{10} h$';
	title_str_Z = sprintf('%s \n %s', line1_Z, line2);
	print_str_Z = 'Z0';
	plot_errs_params_Z = {title_str_Z, print_str_Z};
	[alpha_Z0, beta_Z0] = plot_initializer_errors(stepsizes, Z0_errs, plot_errs_params_Z);

	% U component
	line1_U = '$\log_{10} \| U_{n+1}^{(0)} - U_{n+1} \|$ vs.~$\log_{10} h$';
	title_str_U = sprintf('%s \n %s', line1_U, line2);
	print_str_U = 'U0';
	plot_errs_params_U = {title_str_U, print_str_U};
	[alpha_U0, beta_U0] = plot_initializer_errors(stepsizes, U0_errs, plot_errs_params_U);

else
	error('Error (main): Invalid combination of error_plotting_option and s.')
end




