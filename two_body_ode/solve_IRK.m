function yns = solve_IRK(f, Df, IRK_params)
	% function yns = solve_IRK(f, Df, IRK_params).
	% Sovles the ODE y' = f(y) using the IRK method (b,c,A).  
	%
	% It would be better to use a different implementation if the method is
	% explicit, although this should work for an explicit method.
	%
	% This is meant to work for any s-stage method.  However, this does NOT
	% exploit the special structure of Lobatto IIIA or Lobatto IIIB methods.
	% It instead solves these as any other RK method.
	%
	% Inputs: 
	% 	f: The RHS function of the ODE
	%	Df: Jacobian of f
	%	IRK_params: cell array structure of {h,N,y0,A,b,c,tol}
	%	h: Step size
	%	N: Number of steps to take.
	% 	y0: the starting value
	%	(b,c,A): The RK coefficients.  The vectors b and c need to be 
	%		column vectors of the same dimension s, and A needs to be 
	%		s by s.
	%	tol: tolerance for the iterative method solving the nonlinear
	%		system of equations.
	% 
	% Outputs:
	%	yns: Matrix of steps.  Each column yns(:,i) is the step i.

	% unpack IRK_params
	h = IRK_params{1};
	N = IRK_params{2};
	y0 = IRK_params{3};
	A = IRK_params{4};
	b = IRK_params{5};
	c = IRK_params{6};
	tol = IRK_params{7};

	%% check dimensions and initialize F, DF, yn, yns

	% check that y0, b, and c are column vectors
	[n m] = size(y0);

	if m>1
		error('Error (solve_IRK): y0 is not a column vector.')
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

	yn = y0;
	yns = yn;

	%% main loop

	for i=1:N
		% solve the nonlinear system H(Y)=0 to find the internal stages
		% Y=[Y1;...;Ys].  Let F(Y)=[f(Y1); ... ; f(Ys)].
		% Note that H must be updated at every step because it depends on yn.
		H = @(Y) Y - kron(ones(s,1),yn) - h*kron(A,eye(n))*F(Y);
		DH = @(Y) eye(n*s) - h*kron(A,eye(n))*DF(Y);

		Y0 = kron(ones(s,1),yn); % use a trivial starting approximation
		[Ynp iter] = newton(H,DH,Y0,tol);
		%[Ynp iter] = simplified_newton(H,DH(Y0),Y0,tol);

		ynp = yn + h*kron(b',eye(n))*computeF(Ynp,f,n,s);
		
		yns = [yns ynp];
		yn = ynp;
	end

end

