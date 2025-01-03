function result = find_error_SA(h, params)
	% function result = find_error_SA(h, params)
	% Finds the error in the starting approximations 
	% for the two-body problem.
	%
	% Only two steps are taken: The first with stepsize h,
	% and the second with stepsize rh.
	%
	% Here y = (q,p) \in R^{2d}.
	%
	% Inputs:
	% 	h: stepsize
	% 	params: {d, A, b, c, ynm, Y0_trivial, tol, r, method_str}
	%	d: common dimension of q and p (d=2 for this problem)
	% 	A, b, c: RK coefficients
	%	ynm: This is y_{n-1}, the initial value.
	%	Y0_trivial: initial guess for first step 
	%	tol: tolerance for iterative scheme to solve the nonlinear systems
	%	r: ratio of the step sizes
	%	method_str: String describing the RK method.  Currently the acceptable
	%		values are 'lobattoIIIA', 'lobattoIIIB', and 'gauss'.
	%
	% Outputs: Returns the cell array 
	% 	result = {Ynp0_err, yns, energy_vec}
	%	Ynp0_err: Error in the initializer y_{n+1}^{(0)}
	%	yns: The steps y_n
	%	energy_vec: vector of energy values.

	% unpack params
	d = params{1};
	A = params{2};
	b = params{3};
	c = params{4};
	ynm = params{5};
	Y0_trivial = params{6};
	tol = params{7};
	r = params{8};
	method_str = params{9};

	% keep only lower case method string without spaces
	method_str = lower(method_str); 
	method_str = strrep(method_str, ' ', ''); 

	ny = 2*d; % dimension of y=(q,p)

	%% first step

	% parameters for G and DG
	G_params = {d}; 

	G = @(y) computeG(y, G_params);
	DG = @(y) computeDG(y, G_params);

	[yn Yn] = solve_IRK_single_step(G, DG, h, ynm, Y0_trivial, A, b, c, tol, 'simplified_newton'); 

	%% second step

	SA_params = {h, b, r, G, method_str};
	Ynp0 = find_SA(ynm, Yn, SA_params);

	[ynp Ynp] = solve_IRK_single_step(G, DG, r*h, yn, Ynp0, A, b, c, tol, 'simplified_newton'); 

	%% find the error

	% if the method is lobatto IIIA, disregard the first internal stage error
	if strcmp(method_str, lower('lobattoIIIA')) 
		Ynp0_err = norm(Ynp0(ny+1:end) - Ynp(ny+1:end));
	else
		Ynp0_err = norm(Ynp0-Ynp);
	end

	yns = [ynm yn ynp];

	% make the energy vector
	qns = yns(1:2,:);
	pns = yns(3:4,:);
	for i=1:3
		qn = qns(:,i);
		pn = pns(:,i);
		if i==1
			energy_vec = (1/2)*norm(pn)^2 - 1/norm(qn);
		else
			Enp = (1/2)*norm(pn)^2 - 1/norm(qn);
			energy_vec = [energy_vec Enp];
		end
	end

	result = {Ynp0_err, yns, energy_vec};

end

