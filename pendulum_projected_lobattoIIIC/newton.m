function [root iter] = newton(F,DF,x0,tol)
% function [root iter] = newton(F,DF,x0,tol)
% Solves F(x) = 0 by Newton's method.

sx0 = size(x0);
sF = size(F(x0));
sDF = size(DF(x0));

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
		error('Newtons method did not converge')
	end

	delta = DF(xn) \ (-F(xn));
	xn = xn + delta;
end

if any(isnan(xn))
	error('Error (newton): The root contains one or more NaN values')
else
	root = xn;
end

