function [root iter] = simplified_newton(F,DF0,x0,tol)
% function [root iter] = simplified_newton(F,DF0,x0,tol)
% Solves F(x) = 0 by the simplified Newton's method.
% DF0 = DF(x0).

sx0 = size(x0);
sF = size(F(x0));
sDF = size(DF0);

if sx0(2) > 1
	error('Error (newton): x0 is not a column vector.')
elseif sF(2) > 1
	error('Error (newton): F(x0) is not a column vector.')
elseif sDF(1) ~= sDF(2)
	error('Error (newton): DF(x0) is not a square matrix.')
elseif sDF(2) ~= sF(1)
	error('Error (newton): F(x0) and DF(x0) have incompatible dimensions.')
end

xn = x0;
iter = 0;

while max(abs(F(xn))) > tol
	iter = iter + 1;
	if iter > 10
		error('Simplified Newton did not converge')
	end

	delta = DF0 \ (-F(xn));
	xn = xn + delta;
end

root = xn;

