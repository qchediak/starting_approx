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
	% Also works for Lobatto IIIC (s=3) as long as projection is used.

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

	if strcmp(method, lower('GaussMethod')) && s==3
		beta_star11 = (1/180)*r*(-30 + 3*(-25 + 6*sqrt(15))*r + 2*(-35 + 9*sqrt(15))*r^2);
		beta_star12 = -(1/45)*r*(3*(-5 + sqrt(15)) + 15*(-4 + sqrt(15))*r + (-35 + 9*sqrt(15))*r^2);
		beta_star13 = (1/180)*r*(30*(-4 + sqrt(15)) + 3*(-55 + 14*sqrt(15))*r + 2*(-35 + 9*sqrt(15))*r^2);
		beta_star21 = -(1/72)*r*(6*(5 + sqrt(15)) + 3*(10 + sqrt(15))*r + 10*r^2);
		beta_star22 = (1/18)*r*(6 + 15*r + 5*r^2);
		beta_star23 = (1/72)*r*(6*(-5 + sqrt(15)) + 3*(-10 + sqrt(15))*r - 10*r^2);
		beta_star31 = -(1/180)*r*(30*(4 + sqrt(15)) + 3*(55 + 14*sqrt(15))*r + 2*(35 + 9*sqrt(15))*r^2);
		beta_star32 = (1/45)*r*(3*(5 + sqrt(15)) + 15*(4 + sqrt(15))*r + (35 + 9*sqrt(15))*r^2);
		beta_star33 = -(1/180)*r*(30 + 3*(25 + 6*sqrt(15))*r + 2*(35 + 9*sqrt(15))*r^2);
		beta_star = [beta_star11 beta_star12 beta_star13; 
			beta_star21 beta_star22 beta_star23;
			beta_star31 beta_star32 beta_star33];
		E = [0 0 1; 0 1 0; 1 0 0];
		beta = kron(ones(s,1), b') - beta_star*E;

		% for this method, we could choose beta_hat = beta
		%beta_hat = beta;

		beta_hat_star11 = -(sqrt(15)/3)*r*c(1) - (sqrt(15)/3)*r^2*c(1)^2;
		beta_hat_star12 = (-1 + sqrt(15)/3)*r*c(1) + (sqrt(15)/3)*r^2*c(1)^2;
		beta_hat_star13 = 0;
		beta_hat_star21 = -(sqrt(15)/3)*r*c(2) - (sqrt(15)/3)*r^2*c(2)^2;
		beta_hat_star22 = (-1 + sqrt(15)/3)*r*c(2) + (sqrt(15)/3)*r^2*c(2)^2;
		beta_hat_star23 = 0;
		beta_hat_star31 = -(sqrt(15)/3)*r*c(3) - (sqrt(15)/3)*r^2*c(3)^2;
		beta_hat_star32 = (-1 + sqrt(15)/3)*r*c(3) + (sqrt(15)/3)*r^2*c(3)^2;
		beta_hat_star33 = 0;
		beta_hat_star = [beta_hat_star11 beta_hat_star12 beta_hat_star13; 
			beta_hat_star21 beta_hat_star22 beta_hat_star23; 
			beta_hat_star31 beta_hat_star32 beta_hat_star33];
		beta_hat = kron(ones(s,1), b') - beta_hat_star*E;

		D_star11 = sqrt(15)/3 + (2*sqrt(15)/3)*r*c(1);
		D_star12 = 1-sqrt(15)/3-(2*sqrt(15)/3)*r*c(1);
		D_star13 = 0;
		D_star21 = sqrt(15)/3 + (2*sqrt(15)/3)*r*c(2);
		D_star22 = 1-sqrt(15)/3-(2*sqrt(15)/3)*r*c(2);
		D_star23 = 0;
		D_star31 = sqrt(15)/3 + (2*sqrt(15)/3)*r*c(3);
		D_star32 = 1-sqrt(15)/3-(2*sqrt(15)/3)*r*c(3);
		D_star33 = 0;
		D_star = [D_star11 D_star12 D_star13; 
			D_star21 D_star22 D_star23; 
			D_star31 D_star32 D_star33];
		D = D_star*E;

	elseif strcmp(method, lower('GaussMethod')) && s==2
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
	elseif strcmp(method, lower('LobattoIIIC')) && s==3
		beta_star11 = -4*r^3;
		beta_star12 = 8*r^3;
		beta_star13 = -4*r^3;
		beta_star21 = -r*(r^2+9*r+12);
		beta_star22 = 2*r^2*(r+6);
		beta_star23 = -r^2*(r+3);
		beta_star31 = -4*r*(4*r^2+9*r+6);
		beta_star32 = 16*r^2*(2*r+3);
		beta_star33 = -4*r^2*(4*r+3);
		beta_star = (1/24) * [beta_star11 beta_star12 beta_star13; 
			beta_star21 beta_star22 beta_star23; 
			beta_star31 beta_star32 beta_star33];
		E = [0 0 1;
			0 1 0; 
			1 0 0];
		beta = kron(ones(s,1), b') - beta_star*E;

		theta_star11 = -2*r^2;
		theta_star12 = 4*(r^2-1);
		theta_star13 = 2-2*r^2;
		theta_star21 = r^2;
		theta_star22 = -5*r^2 - 12*r - 4;
		theta_star23 = 4*r^2 + 6*r + 2;
		theta_star31 = -2*r^2;
		theta_star32 = -4*(2*r^2 + 6*r + 1);
		theta_star33 = 2*(5*r^2 + 6*r + 1);

		theta_star = (1/12) * [theta_star11 theta_star12 theta_star13;
								theta_star21 theta_star22 theta_star23;
								theta_star31 theta_star32 theta_star33];

		e1 = [1; 0; 0];
		beta_hat_star = theta_star + b(s)*kron(ones(s,1),e1');
		beta_hat = kron(ones(s,1), b') - beta_hat_star*E;

		D_star11 = 0;
		D_star12 = 2;
		D_star13 = -1;

		D_star21 = 0;
		D_star22 = r+2;
		D_star23 = -r-1;

		D_star31 = r;
		D_star32 = 2;
		D_star33 = -r-1;

		D_star = [D_star11 D_star12 D_star13; 
			D_star21 D_star22 D_star23;
			D_star31 D_star32 D_star33];

		D = D_star*E;
	else
		error('Error (find_SA_DAE): Invalid method and s value combination')
	end

	%% compute Ynp0, Znp0, Unp0

	Ynp0 = kron(ones(s,1), ynm) + h * kron(beta, eye(ny)) * computeF(Yn,Zn,f,ny,s);
	Znp0 = kron(ones(s,1), znm) + h * kron(beta_hat, eye(nz)) * computeK(Yn,Zn,Un,k,nz,nu,s);
	Unp0 = kron(D, eye(nu)) * Un;
end


