function [ynp Ynp] = solve_IRK_single_step(f, Df, h, yn, Y0, A, b, c, tol, iterative_scheme)
	%function ynp = solve_IRK_single_step(f, Df, h, yn, y_SA, A, b, c, tol, iterative_scheme)
	% Sovles the ODE y' = f(y) using the IRK method (b,c,A) over a
	% single step.  
	%
	% It would be better to use a different implementation if the method is
	% explicit, although this should work for an explicit method.
	%
	% This is meant to work for any s-stage method.  However, this does NOT
	% exploit the special structure of Lobatto IIIA or Lobatto IIIB methods.
	% It instead solves these as it would any other IRK method.
	%
	% Inputs: 
	% 	f: The RHS function of the ODE
	%	Df: Jacobian of f
	%	h: Step size
	% 	yn: the starting value
	%	Y0: The starting approximation used in solving the nonlinear
	%		system
	%	(b,c,A): The RK coefficients.  The vectors b and c need to be 
	%		column vectors of the same dimension s, and A needs to be 
	%		s by s.
	%	tol: tolerance for the iterative method solving the nonlinear
	%		system of equations.
	%	iterative_scheme: The iterative scheme used to solve the
	%		nonlinear system.  Can take either of the values 'newton' or
	%		'simplified_newton'.
	% 
	% Outputs:
	%	ynp: The updated step.
	%	Ynp: Vector of internal stages. 

	%% check dimensions and initialize F, DF, yn, yns

	% check that yn, b, and c are column vectors
	[n m] = size(yn);

	if m>1
		error('Error (solve_IRK): yn is not a column vector.')
	end

	% check that the dimensions of b,c,A are compatible
	[s m] = size(b);
	if m>1
		error('Error (solve_IRK): b is not a column vector.')
	elseif size(c,1) ~= s
		error('Error (solve_IRK): c is not the same length as b')
	elseif size(c,2) ~= 1
		error('Error (solve_IRK): c is not a column vector')
	elseif size(A,1) ~= s || size(A,2) ~= s
		error('Error (solve_IRK): A is not s by s')
	end

	% for more readable code in the main loop
	F = @(Y) computeF(Y,f,n,s);
	DF = @(Y) computeDF(Y,Df,n,s);

	%% main loop

	% solve the nonlinear system H(Y)=0 to find the internal stages
	% Y=[Y1;...;Ys].  Let F(Y)=[f(Y1); ... ; f(Ys)].
	% Note that H must be updated at every step because it depends on yn.
	H = @(Y) Y - kron(ones(s,1),yn) - h*kron(A,eye(n))*F(Y);
	DH = @(Y) eye(n*s) - h*kron(A,eye(n))*DF(Y);

	if strcmp(iterative_scheme, 'newton')
		[Ynp iter] = newton(H,DH,Y0,tol);
	elseif strcmp(iterative_scheme, 'simplified_newton')
		[Ynp iter] = simplified_newton(H,DH(Y0),Y0,tol);
	else
		error('Error (solve_IRK_single_step): Valid values for iterative_scheme are newton and simplified_newton')
	end

	ynp = yn + h*kron(b',eye(n))*F(Ynp);

end

