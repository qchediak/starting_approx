function Xnp0 = find_SA(xnm, Xn, SA_params)
	% function Xnp0 = find_SA(xnm, Xn, SA_params)
	% Finds the starting approximations Xnp0 to Xnp based on the previous 
	% internal stages Xn and the past step xnm.  Assumes we are solving 
	% the ODE x'=f(x).
	%
	% Inputs:
	%	xnm: The step x_{n-1}
	%	Xn: vector of the internal stages X_{n,1},...,X_{n,s}.
	%	SA_params = {h, b, r, f, method}
	%	h: stepsize of the first step (does NOT factor in r)
	%	b: weights of the RK method
	%	r: ratio of second stepsize / first stepsize
	% 	f: RHS of the ODE x'=f(x).
	%	method: String describing the method.  Currently the acceptable
	%		values are 'lobattoIIIA', 'lobattoIIIB', or 'gauss method'.
	%
	% Outputs:
	% 	Xnp0 = [X_{n+1,1}^(0); ...; X_{n+1,s}^(0)]
	%
	% Note params has a cell array structure.
	%
	% Currently only works with Lobatto IIIA, Lobatto IIIB, or Gauss, all
	% with s=2.

	% unpack params
	h = SA_params{1};
	b = SA_params{2};
	r = SA_params{3};
	f = SA_params{4};
	method = SA_params{5};

	nx = size(xnm,1);
	s = length(b);

	%% find starred coefficients

	method = lower(method); % lower case
	method = strrep(method, ' ', ''); % remove white spaces

	if strcmp(method, lower('lobattoIIIA')) && s==2
		beta_star = (1/2) * [0 0; -2*r-r^2 r^2];
		%beta = (1/2) * [1 1; 1-r^2 1+2r+r^2];
		E = [0 1; 1 0];
	elseif strcmp(method, lower('lobattoIIIB')) && s==2
		beta_star = [-(1/2)*r 0; -(1/2)*r 0];
		E = [0 1; 1 0];
	elseif strcmp(method, lower('GaussMethod')) && s==2
		beta_star11 = -(1/12)*r*(2*sqrt(3) + (-3+2*sqrt(3))*r);
		beta_star12 = (1/12)*(-3+2*sqrt(3))*r*(2+r);
		beta_star21 = -(1/12)*(3+2*sqrt(3))*r*(2+r);
		beta_star22 = (1/12)*r*(2*sqrt(3)+(3+2*sqrt(3))*r);
		beta_star = [beta_star11 beta_star12; beta_star21 beta_star22];
		E = [0 1; 1 0];
	else
		error('Error (find_SA): Invalid method and s value combination')
	end

	%% find unstarred coefficients

	beta = kron(ones(s,1), b') - beta_star*E;

	%% compute Xnp0

	F = @(X) computeF(X,f,nx,s); % F(Xn) = [f(X_{n,1});...;f(X_{n,s})]
	Xnp0 = kron(ones(s,1), xnm) + h * kron(beta, eye(nx)) * F(Xn);

end


