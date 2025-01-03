function result = find_error_SA(h, params)
% function result = find_error_SA(h, params)
% Finds the error in the starting approximations 
% for the pendulum problem.
%
% Only two steps are taken: The first with stepsize h,
% and the second with stepsize rh.
%
% Inputs:
% 	h: stepsize
%	params: {ny, nz, nu, yn, zn, g0, l, A, Ahat, b, c, N, U10, U20, U30, tol, r}
%	ny, nz, nu: dimensions of y,z,u
% 	yn, zn: step
%	g0: gravitational constant
% 	l: length of rod
% 	A, Ahat, b, c: PRK coefficients
% 	N: number of steps to take
%	U10, U20, U30: initial guess for first step for u-component
%	tol: tolerance for iterative scheme to solve the nonlinear systems
%
% Outputs: Returns the cell array 
% 	result = {Y0_err, Z0_err, znp0_err, U10_err, U20_err, U30_err, yns, zns, energy_vec};
%	Y0_err: Error in the initializer [Y_{n+1,2}^{(0)}; Y_{n+1,3}^{(0)}]
%	Z0_err: Error in the initializer [Z_{n+1,1}^{(0)}; Z_{n+1,2}^{(0)}; Z_{n+1,3}^{(0)}]
%	znp0_err: Error in the initializer z_{n+1}^{(0)}
%	U10_err: Error in the initializer U_{n+1,1}^{(0)}
%	U20_err: Error in the initializer U_{n+1,2}^{(0)}
%	U30_err: Error in the initializer U_{n+1,3}^{(0)}
%	yns, zns: matrices of steps.
%	energy_vec: vector of energy values.

% unpack params
ny = params{1};
nz = params{2};
nu = params{3};
yn = params{4};
zn = params{5};
g0 = params{6};
l = params{7};
A = params{8};
Ahat = params{9};
b = params{10};
c = params{11};
U10_trivial = params{12};
U20_trivial = params{13};
U30_trivial = params{14};
tol = params{15};
r = params{16};

% matrices of yn and zn values
yns = yn;
zns = zn;

% Initial energy, where potential energy is zero when the pendulum is at
% its lowest point.  Initially (yn,zn) = (y0,z0).
E = (1/2)*norm(zn)^2 + g0*(l+yn(2)); 

for n=1:2
	% initial guesses for the simplified Newton method.
	if n==1
		x0H = [yn; yn; zn; zn; zn; U10_trivial; U20_trivial];
		x0T = [zn; U30_trivial];
	else
		Yn1 = ynm;
		Yn2 = Y2;
		Yn3 = Y3;
		Zn1 = Z1;
		Zn2 = Z2;
		Zn3 = Z3;
		Un1 = U1;
		Un2 = U2;
		Un3 = U3;

		xnm = {ynm,znm};
		Xn = {Yn1,Yn2,Yn3,Zn1,Zn2,Zn3,Un1,Un2,Un3};
		SA_params = {ny,nz,nu,h,yn,zn,g0,l,A,Ahat,b,c,r};

		Xnp0 = find_SA(xnm, Xn, SA_params);

		Ynp20 = Xnp0{1};
		Ynp30 = Xnp0{2};
		Znp10 = Xnp0{3};
		Znp20 = Xnp0{4};
		Znp30 = Xnp0{5};
		znp0 = Xnp0{6};
		Unp10 = Xnp0{7};
		Unp20 = Xnp0{8};
		Unp30 = Xnp0{9};

		x0H = [Ynp20; Ynp30; Znp10; Znp20; Znp30; Unp10; Unp20];
		x0T = [znp0; Unp30];
	end

%%	%% solve the first nonlinear system
	% The function H needs to be redefined at every step using Fparams.
	% The argument of H is x = [Ybar; Z; Utilde], where
	% Ybar = [Ynp2; Ynp3], Z = [Znp1; Znp2; Znp3], Utilde = [Unp1; Unp2].

	% Step size is h for the first step and rh for the second.
	% Note this changes Gparams as well.

	A0 = A(2:3,:);
	A0hat = Ahat(:,1:2);

	if n==1
		Hparams = {ny, nz, nu, h, yn, zn, g0, l, A0, A0hat};
	else
		Hparams = {ny, nz, nu, r*h, yn, zn, g0, l, A0, A0hat};
	end

	H = @(x) H1(x,Hparams);
	DH = @(x) DH1(x,Hparams);
	[H_root H_iter] = simplified_newton(H,DH(x0H),x0H,tol);
	%[H_root H_iter] = newton(H,DH,x0H,tol);

	% unpack H_root
	Ybar = H_root(1:2*ny);
	Z = H_root(2*ny+1:2*ny+3*nz);
	Utilde = H_root(2*ny+3*nz+1:2*ny+3*nz+2*nu);

	% unpack Ybar
	Y2 = Ybar(1:ny);
	Y3 = Ybar(ny+1:2*ny);

	% unpack Z
	Z1 = Z(1:nz);
	Z2 = Z(nz+1:2*nz);
	Z3 = Z(2*nz+1:3*nz);

	% unpack Utilde
	U1 = Utilde(1:nu)
	U2 = Utilde(nu+1:2*nu)

%%	%% solve the second nonlinear system
	% The function T needs to be redefined at every step using Tparams.

	% Note that Tparams is Tparams with three new parameters.  This
	% makes the code easier to read and write but introduces some
	% inefficiency because not all of these parameters are actually
	% used.

	% The argument of T is x = [znp; Unp3].

	Tparams = [Hparams {b, Ybar, Utilde}];
	T = @(x) T1(x,Tparams);
	DT = @(x) DT1(x,Tparams);
	[T_root T_iter] = simplified_newton(T,DT(x0T),x0T,tol);
	%[T_root T_iter] = newton(T,DT,x0T,tol);

	% unpack T_root 
	znp = T_root(1:nz);
	U3 = T_root(nz+1:nz+nu)

	% Update y
	ynp = Y3;

	% update the matrix of yns, zns
	yns = [yns ynp];
	zns = [zns znp];

	% update energy
	Enp = (1/2)*norm(znp)^2 + g0*(l+ynp(2)); 
	E = [E Enp];

	ynm = yn;
	znm = zn;

	yn = ynp;
	zn = znp;
end	

energy_vec = E; % to avoid confusion with error

%% find error

% find error in the starting approximations used in the second step
Y0_err = norm([Ynp20; Ynp30] - [Y2; Y3]);
Z0_err = norm([Znp10; Znp20; Znp30] - [Z1; Z2; Z3]);
znp0_err = norm(znp0 - znp);
U10_err = norm(Unp10 - U1);
U20_err = norm(Unp20 - U2);
U30_err = norm(Unp30 - U3);

result = {Y0_err, Z0_err, znp0_err, U10_err, U20_err, U30_err, yns, zns, energy_vec};

