function [znp Unpsp iter] = solve_second_system(params, h, zn, Ynp, Znp, Unp, x0, iterative_scheme)
	% function [znp Unpsp iter2] = solve_second_system(params, h, zn, Ynp, Znp, Unp, x0, iterative_scheme)
	%
	% Solves the second nonlinear system in the projected Lobatto IIIC
	% method.
	%
	% Inputs:
	%	params: Cell array consisting of: g0, l, tol, ny, nz, nu,
	%		method_str, b, c, A, Ahat, f, k, g, gyf.
	%	h: step size.
	%	zn: z-component of the update step from the most recent step.
	%	Ynp, Znp, Unp: Vectors of internal stages already found in the
	%		first nonlinear system.
	%	x0 = (znp0, Unpsp0): The starting approximation for the iterative
	%		solver for the second system of nonlinear equations.
	%	iterative_scheme: The iterative scheme used to solve the
	%	nonlinear system.  Can take the following values: 
	%		'newton': The standard Newton method.
	%		'simplified newton': Uses the matrix DF(x0)
	%
	% Outputs:
	%	znp: The z-component of the update step
	% 	Unpsp: The stage U_{n+1,s+1}
	% 	iter: Number of iterations used by the iterative scheme.

	% unpack params
	g0 = params{1};
	l = params{2};
	tol = params{3};
	ny = params{4};
	nz = params{5};
	nu = params{6};
	method_str = params{7};
	b = params{8};
	c = params{9};
	A = params{10};
	Ahat = params{11};
	f = params{12};
	k = params{13};
	g = params{14};
	gyf = params{15};

	s = size(A,1);

	% ---------------------------
	% Solve the nonlinear system
	% ---------------------------

	% The function T needs to be redefined at every step using Tparams.
	% The argument of T is x = [znp; Unpsp].

	Tparams = {ny, nz, nu, h, zn, Ynp, Znp, Unp, b, k, gyf};

	if strcmpi(iterative_scheme, 'newton')
		T = @(x) T1(x,Tparams);
		DT = @(x) DT1(x,Tparams);
		[root iter] = newton(T,DT,x0,tol);
	elseif strcmpi(iterative_scheme, 'simplified newton')
		T = @(x) T1(x,Tparams);
		DT = @(x) DT1(x,Tparams);
		[root iter] = simplified_newton(T,DT(x0),x0,tol);
	else
		error('Error (solve_second_system): Unknown iterative scheme');
	end

	% unpack root
	znp = root(1:nz);
	Unpsp = root(nz+1:nz+nu);
end

