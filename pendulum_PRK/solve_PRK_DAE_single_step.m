function [ynp znp Ynp Znp Unp iter] = solve_PRK_DAE_single_step(params, h, yn, zn, un, x0, iterative_scheme)
	% function [ynp znp Ynp Znp Unp] = solve_IRK_DAE_single_step(params, h, yn, zn, x0, iterative_scheme)
	%
	% Solves the index 3 DAE pendulum problem: y'=f(y,z), z'=k(y,z,u), 0=g(y).
	%
	% This is specific to the pendulum problem; however, it could be
	% modified for other index 3 DAEs.
	%
	% This is meant to work with any s-stage PRK method.  However, it
	% does NOT exploit the special structure of the Lobatto IIIA-IIIB
	% method, so using this function with Lobatto IIIA-IIIB is
	% suboptimal.
	%
	% Inputs: 
	%	params: Cell array consisting of: g0, l, tol, ny, nz, nu,
	%		method_str, b, c, A, Ahat, f, k, g
	%	h: step size.
	%	(yn, zn, un): The most recent step.  If desired, un could be taken
	%		to be U_s.
	%	x0 = (Y0,Z0,U0): The starting approximation for the iterative
	%		solver for the system of nonlinear equations.
	%	iterative_scheme: The iterative scheme used to solve the
	%	nonlinear system.  Can take the following values: 
	%		'newton': The standard Newton method.
	%		'simplified newton 1': Uses the matrix DF(x0)
	%		'simplified newton 2': Uses a matrix J which exploits the
	%			structure of the index 3 problem and approximates h as small.
	%
	% Outputs:
	%	ynp, znp: The update step
	% 	Ynp, Znp, Unp: Vectors of the internal stages
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

	s = size(A,1);

	% ---------------------------
	% Solve the nonlinear system
	% ---------------------------

	% The function H needs to be redefined at every step using Hparams.
	% The argument of H is x = [Y; Z; U],
	% where Y = [Y1;...;Ys], Z = [Z1;...;Zs], U = [U1;...;Us].

	Hparams = {ny, nz, nu, h, yn, zn, g0, l, A, Ahat, f, k, g};

	if strcmpi(iterative_scheme, 'newton')
		H = @(x) H1(x,Hparams);
		DH = @(x) DH1(x,Hparams);
		[root iter] = newton(H,DH,x0,tol);
	elseif strcmpi(iterative_scheme, 'simplified newton 1')
		H = @(x) H1(x,Hparams);
		DH = @(x) DH1(x,Hparams);
		[root iter] = simplified_newton(H,DH(x0),x0,tol);
	elseif strcmpi(iterative_scheme, 'simplified newton 2')
		H = @(x) H2(x,Hparams); % alternate formulation of the nonlinear system
		J = computeJ(un,Hparams); % can be used in simplified Newton
		[root iter] = simplified_newton(H,J,x0,tol);
	else
		error('Error (solve_IRK_DAE_single_step): Unknown iterative scheme');
	end

	% unpack H_root
	Ynp = root(1:ny*s);			% for s=2, ny=nz=2, nu=1, this is 1:4
	Znp = root(ny*s+1:ny*s+nz*s);	% for s=2, ny=nz=2, nu=1, this is 5:8
	Unp = root(ny*s+nz*s+1:end);	% for s=2, ny=nz=2, nu=1, this is 9:10

	% Find the update step
	% Note that this assumes ny=nz
	ynp = yn + h*kron(b',eye(ny))*computeF(Ynp,Znp,f,ny,s);
	znp = zn + h*kron(b',eye(nz))*computeK(Ynp,Znp,Unp,k,ny,nu,s);
end

