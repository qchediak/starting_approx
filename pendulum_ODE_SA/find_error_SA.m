function result = find_error_SA(h, params)
	% function result = find_error_SA(h, params)
	% Finds the error in the starting approximations 
	% for the pendulum problem.
	%
	% Only two steps are taken: The first with stepsize h,
	% and the second with stepsize rh.
	%
	% Here psi = (y,z,u).
	%
	% Inputs:
	% 	h: stepsize
	%	params: {ny, nz, nu, g0, l, A, b, c, psinm, psi0_trivial, tol, r, method}
	%	ny, nz, nu: dimensions of y,z,u
	%	g0: gravitational constant
	% 	A, b, c: RK coefficients
	%	psinm: This is psi_{n-1}, the initial value.
	%	Psi0_trivial: initial guess for first step 
	%	tol: tolerance for iterative scheme to solve the nonlinear systems
	%	r: ratio of the step sizes
	%	method: String describing the RK method.  Currently the acceptable
	%		values are 'lobattoIIIA' and 'lobattoIIIB'.
	%
	% Outputs: Returns the cell array 
	% 	result = {Psinp0_err, psins, energy_vec}
	%	Psinp0_err: Error in the initializer psi_{n+1}^{(0)}
	%	psins: The steps psi_n
	%	energy_vec: vector of energy values.

	% unpack params
	ny = params{1};
	nz = params{2};
	nu = params{3};
	g0 = params{4};
	l = params{5};
	A = params{6};
	b = params{7};
	c = params{8};
	psinm = params{9};
	Psi0_trivial = params{10};
	tol = params{11};
	r = params{12};
	method = params{13};

	% keep only lower case method string without spaces
	method = lower(method); 
	method = strrep(method, ' ', ''); 

	npsi = ny + nz + nu; % dimension of psi=(y,z,u)

	%% first step

	% parameters for G and DG
	G_params = {g0, ny, nz, nu}; 

	G = @(psi) computeG(psi, G_params);
	DG = @(psi) computeDG(psi, G_params);

	[psin Psin] = solve_IRK_single_step(G, DG, h, psinm, Psi0_trivial, A, b, c, tol, 'simplified_newton'); 

	%% second step

	SA_params = {h, b, r, G, method};
	Psinp0 = find_SA(psinm, Psin, SA_params);

	[psinp Psinp] = solve_IRK_single_step(G, DG, r*h, psin, Psinp0, A, b, c, tol, 'simplified_newton'); 

	%% find the error

	% if the method is lobatto IIIA, disregard the first internal stage error
	if strcmp(method, lower('lobattoIIIA')) 
		Psinp0_err = norm(Psinp0(npsi+1:end) - Psinp(npsi+1:end));
	else
		Psinp0_err = norm(Psinp0-Psinp);
	end

	psins = [psinm psin psinp];

	% make the energy vector
	yns = psins(1:2,:);
	zns = psins(3:4,:);
	for i=1:3
		yn = yns(:,i);
		zn = zns(:,i);
		if i==1
			energy_vec = (1/2)*norm(zn)^2 + g0*(l+yn(2));
		else
			Enp = (1/2)*norm(zn)^2 + g0*(l+yn(2));
			energy_vec = [energy_vec Enp];
		end
	end

	result = {Psinp0_err, psins, energy_vec};

end

