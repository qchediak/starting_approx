function [znp0 Unpsp0] = find_SA_DAE_second_system(ynm, znm, Yn, Zn, Un, SA_params)
	% function [znp0 Unpsp0] = find_SA_DAE_second_system(ynm, znm, Yn, Zn, Un, SA_params)
	%
	% Finds starting approximations for the second nonlinear system of
	% the projected lobatto IIIC method applied to an index 3 DAE.  For
	% simplicity, the inputs are chosen to match those of the function
	% find_SA_DAE; they are not all actually used.
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
	%	method: String describing the method.  
	%
	% Outputs:
	% 	znp0 = z_{n+1}^{(0)}
	% 	Unpsp0 = U_{n+1,s+1}^{(0)}
	%
	% Note SA_params has a cell array structure.

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

	method = lower(method); % lower case
	method = strrep(method, ' ', ''); % remove white spaces

	% Check that method and s are as expected
	if strcmp(method, lower('LobattoIIIC')) && s==3
		e1 = eye(1,s);
		E = [0 0 1; 0 1 0; 1 0 0];

		theta_star41 = 0;
		theta_star42 = -r^2 - 2*r - 1/3;
		theta_star43 = r^2 + r + 1/6;

		theta_star4 = [theta_star41 theta_star42 theta_star43];
		beta_hat_star4 = theta_star4 + b(s)*e1;
		beta_hat4 = b' - beta_hat_star4*E;

		D_star41 = 0;
		D_star42 = 2*r+2;
		D_star43 = -2*r-1;

		D_star4 = [D_star41 D_star42 D_star43];
		D4 = D_star4*E;
	else
		error('Error (find_SA_DAE_second_system): Invalid method and s value combination')
	end

	% Compute znp0 and Unpsp0

	znp0 = znm + h * kron(beta_hat4, eye(nz)) * computeK(Yn,Zn,Un,k,nz,nu,s);
	Unpsp0 = kron(D4, eye(nu)) * Un;
end


