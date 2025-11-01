function [root iter] = simplified_newton(F,J,x0,tol)
	% function [root iter] = simplified_newton(F,J,x0,tol)
	%
	% Solves F(x) = 0 by the simplified Newton's method, using a
	% constant matrix J which approximates the Jacobian.  Can be taken
	% to be J=DF(x0).
	%
	% Inputs:
	%	F: The nonlinear system to solve is F(x)=0.
	%	J: The approximate Jacobian used in the step delta = J \ (-F(xn))
	%	x0: The initial guess
	%	tol: Tolerance.
	%
	% Outputs:
	%	root: The root of F(x)=0
	%	iter: The number of iterations to find the root.

	if size(x0,2) > 1
		error('Error (simplified_newton): x0 is not a column vector.')
	elseif size(F(x0),2) > 1
		error('Error (simplified_newton): F(x0) is not a column vector.')
	elseif size(J,1) ~= size(J,2)
		error('Error (simplified_newton): J is not a square matrix.')
	elseif size(J,2) ~= size(F(x0),1)
		error('Error (simplified_newton): F(x0) and J have incompatible dimensions.')
	end

	xn = x0;
	iter = 0;

	while max(abs(F(xn))) > tol
		iter = iter + 1;
		if iter > 10
			error('simplified_newton did not converge')
		end

		delta = J \ (-F(xn));
		xn = xn + delta;
	end

	if any(isnan(xn))
		error('Error (simplified_newton): The root contains one or more NaN values')
	else
		root = xn;
	end
end

