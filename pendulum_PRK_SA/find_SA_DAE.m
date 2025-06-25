function [Ynp0 Znp0 Unp0] = find_SA_DAE(ynm, znm, Yn, Zn, Un, SA_params)
	% function [Ynp0 Znp0 Unp0] = find_SA_DAE(ynm, znm, Yn, Zn, Un, SA_params)
	%
	% Finds the starting approximations Ynp0, Znp0, Unp0 to Ynp, Znp, Unp based
	% on the previous internal stages Yn, Zn, Un and the past step ynm, znm.  	
	%
	% Inputs:
	%	ynm, znm: The past step y_{n-1}, z_{n-1}
	%	Yn: vector of the internal stages Y_{n,1},...,Y_{n,s}.
	%	Zn: vector of the internal stages Z_{n,1},...,Z_{n,s}.
	%	Un: vector of the internal stages U_{n,1},...,U_{n,s}.
	%	SA_params = {h, b, c, r, f, k, method}
	%	h: stepsize of the first step (does NOT factor in r)
	%	b: weights of the PRK method
	%	c: nodes of the PRK method
	%	r: ratio of second stepsize / first stepsize
	% 	f: RHS of y'=f(y,z).
	% 	k: RHS of z'=f(y,z,u).
	%	nu: dimension of u
	%	method: String describing the method.  Currently the acceptable
	%		values are 'Gauss', or 'Radau IIA' (spacing and
	%		capitalization don't matter)
	%
	% Outputs:
	% 	Ynp0 = [Y_{n+1,1}^(0); ...; Y_{n+1,s}^(0)]
	% 	Znp0 = [Z_{n+1,1}^(0); ...; Z_{n+1,s}^(0)]
	% 	Unp0 = [U_{n+1,1}^(0); ...; U_{n+1,s}^(0)]
	%
	% Note params has a cell array structure.
	%
	% Currently works with Gauss or Radau IIA, both with s=2.

	% unpack params
	h = SA_params{1};
	b = SA_params{2};
	c = SA_params{3};
	r = SA_params{4};
	f = SA_params{5};
	k = SA_params{6};
	nu = SA_params{7};
	method = SA_params{8};

	ny = size(ynm,1);
	nz = size(znm,1);
	s = length(b);

	%% find starred coefficients

	method = lower(method); % lower case
	method = strrep(method, ' ', ''); % remove white spaces

	if strcmp(method, lower('Gauss')) && s==2
		beta_star11 = -(1/12)*r*(2*sqrt(3) + (-3+2*sqrt(3))*r);
		beta_star12 = (1/12)*(-3+2*sqrt(3))*r*(2+r);
		beta_star21 = -(1/12)*(3+2*sqrt(3))*r*(2+r);
		beta_star22 = (1/12)*r*(2*sqrt(3)+(3+2*sqrt(3))*r);
		beta_star = [beta_star11 beta_star12; beta_star21 beta_star22];
		E = [0 1; 1 0];
		beta = kron(ones(s,1), b') - beta_star*E;

		beta_hat11 = 1+r*c(1);
		beta_hat12 = 0;
		beta_hat21 = 1+r*c(2);
		beta_hat22 = 0;
		beta_hat = [beta_hat11 beta_hat12; beta_hat21 beta_hat22];

		D11 = 1;
		D12 = 0;
		D21 = 1;
		D22 = 0;
		D = [D11 D12; D21 D22];
	elseif strcmp(method, lower('RadauIIA')) && s==2
		beta_star11 = (1/12) * (-r) * (r+4);
		beta_star12 = (1/12) * r^2;
		beta_star21 = (1/12) * (-3*r) * (3*r+4);
		beta_star22 = (1/12) * 9 * r^2;
		beta_star = [beta_star11 beta_star12; beta_star21 beta_star22];
		E = [0 1; 1 0];
		beta = kron(ones(s,1), b') - beta_star*E;

		beta_hat11 = 1+r*c(1);
		beta_hat12 = 0;
		beta_hat21 = 1+r*c(2);
		beta_hat22 = 0;
		beta_hat = [beta_hat11 beta_hat12; beta_hat21 beta_hat22];

		D11 = 1;
		D12 = 0;
		D21 = 1;
		D22 = 0;
		D = [D11 D12; D21 D22];
	else
		error('Error (find_SA): Invalid method and s value combination')
	end

	%% compute Ynp0, Znp0, Unp0

	Ynp0 = kron(ones(s,1), ynm) + h * kron(beta, eye(ny)) * computeF(Yn,Zn,f,ny,s);
	Znp0 = kron(ones(s,1), znm) + h * kron(beta_hat, eye(nz)) * computeK(Yn,Zn,Un,k,nz,nu,s);
	Unp0 = kron(D, eye(nu)) * Un;
end


