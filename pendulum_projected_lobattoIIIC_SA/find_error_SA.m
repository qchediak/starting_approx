function result = find_error_SA(h, params)
	% function result = find_error_SA(h, params)
	% Finds the error in the starting approximations 
	% for the pendulum problem.
	%
	% Only two steps are taken: The first with stepsize h,
	% and the second with stepsize rh.
	%
	% Inputs:
	% 	h: stepsize
	%	params: Cell array of ny, nz, nu, g0, l, b, c, A, Ahat, ynm, znm, unm,
	%		x0_trivial, tol, r, method, f, k, g, gyf, iterative_scheme1,
	%		iterative_scheme2
	%	ny, nz, nu: dimensions of y,z,u
	%	g0: gravitational constant
	% 	b, c, A, Ahat: PRK coefficients
	%	ynm, znm: The initial values (y_{n-1}, z_{n-1})
	%	x0_trivial: initial guess for first step 
	%	tol: tolerance for iterative scheme to solve the nonlinear systems
	%	r: ratio of the step sizes
	%	method: String describing the PRK method.  
	%	f, k, g: The DAE is y'=f(y,z), z'=k(y,z,u), 0=g(y)
	%	gyf: The hidden constraint function g_y(y)f(y,z)
	%	iterative_scheme1, iterative_scheme2: Strings which indicate
	%		which iterative schemes are to be used in solving the two
	%		nonlinear systems.
	%
	% Outputs: Returns the cell array 
	% 	result = {Y0_err_vec, Z0_err_vec, U0_err_vec, znp0_err_vec, Unpsp0_err_vec, yns, zns, energy_vec}
	%	Ynp0_err_vec: Vector of errors in the initializer Y_{n+1}^{(0)}
	%	Znp0_err_vec: Vector of errors in the initializer Z_{n+1}^{(0)}
	%	Unp0_err_vec: Vector of errors in the initializer U_{n+1}^{(0)}
	%	znp0_err_vec: Vector of errors in the initializer z_{n+1}^{(0)}
	%	Unpsp0_err_vec: Vector of errors in the initializer U_{n+1,s+1}^{(0)}
	%	yns: The steps [y_{n-1}, y_n, y_{n+1}]
	%	zns: The steps [z_{n-1}, z_n, z_{n+1}]
	%	energy_vec: vector of energy values.

	% unpack params
	ny = params{1};
	nz = params{2};
	nu = params{3};
	g0 = params{4};
	l = params{5};
	b = params{6};
	c = params{7};
	A = params{8};
	Ahat = params{9};
	ynm = params{10};
	znm = params{11};
	unm = params{12};
	x0_trivial = params{13};
	tol = params{14};
	r = params{15};
	method = params{16};
	f = params{17};
	k = params{18};
	g = params{19};
	gyf = params{20};
	iterative_scheme1 = params{21};
	iterative_scheme2 = params{22};

	% keep only lower case method string without spaces
	method = lower(method); 
	method = strrep(method, ' ', ''); 

	alpha0 = [0; g0];

	% ------------
	% first step
	% ------------

	% parameters for solve_PRK_DAE_single_step
	solve_params = {g0, l, tol, ny, nz, nu, method, b, c, A, Ahat, f, k, g, gyf};

	% first_nonlinear_system
	[yn zn Yn Zn Un iter1] = solve_PRK_DAE_single_step(solve_params, h, ynm, znm, unm, x0_trivial, iterative_scheme1);

	% Second nonlinear system
	[zn Unsp iter2] = solve_second_system(solve_params, h, znm, Yn, Zn, Un, [znm; unm], iterative_scheme2);

	% -----------
	% second step
	% -----------

	un = (1/l^2) * (dot(zn,zn) - dot(yn,alpha0)); % satisfies hidden constraint.  

	SA_params = {h, b, c, r, f, k, nu, method};
	[Ynp0 Znp0 Unp0] = find_SA_DAE(ynm, znm, Yn, Zn, Un, SA_params);
	[znp0 Unpsp0] = find_SA_DAE_second_system(ynm, znm, Yn, Zn, Un, SA_params);

	% First nonlinear system
	[ynp znp Ynp Znp Unp iter3] = solve_PRK_DAE_single_step(solve_params, r*h, yn, zn, un, [Ynp0; Znp0; Unp0], iterative_scheme1);

	% Second nonlinear system
	[znp Unpsp iter4] = solve_second_system(solve_params, r*h, zn, Ynp, Znp, Unp, [znp0; Unpsp0], iterative_scheme2);

	%% find the error

	Ynp0_err_vec = abs(Ynp0 - Ynp);
	Znp0_err_vec = abs(Znp0 - Znp);
	Unp0_err_vec = abs(Unp0 - Unp);

	znp0_err_vec = abs(znp0 - znp);
	Unpsp0_err_vec = abs(Unpsp0 - Unpsp);

	yns = [ynm yn ynp];
	zns = [znm zn znp];

	% make the energy vector
	energy_vec = (1/2)*norm(znm)^2 + g0*(l+ynm(2));
	for i=2:3
		y_vec = yns(:,i);
		z_vec = zns(:,i);
		Enp = (1/2)*norm(z_vec)^2 + g0*(l+y_vec(2));
		energy_vec = [energy_vec Enp];
	end

	result = {Ynp0_err_vec, Znp0_err_vec, Unp0_err_vec, znp0_err_vec, Unpsp0_err_vec, yns, zns, energy_vec};
end

