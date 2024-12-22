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
%	params: {ny, nz, nu, yn, zn, g0, l, A, Ahat, b, c, N, U10, U20, tol, r}
%	ny, nz, nu: dimensions of y,z,u
% 	yn, zn: step
%	g0: gravitational constant
% 	l: length of rod
% 	A, Ahat, b, c: PRK coefficients
% 	N: number of steps to take
%	U10, U20: initial guess for first step for u-component
%	tol: tolerance for iterative scheme to solve the nonlinear systems
%
% Outputs: Returns the cell array 
% 	result = {Y0_err, Z0_err, znp0_err, U0_err, yns, zns, energy_vec};
%	Y0_err: Error in the initializer Y_{n+1,2}^{(0)}
%	Z0_err: Error in the initializer [Z_{n+1,1}^{(0)}; Z_{n+1,2}^{(0)}]
%	znp0_err: Error in the initializer z_{n+1}^{(0)}
%	U0_err: Error in the initializer [U_{n+1,1}^{(0)}; U_{n+1,2}^{(0)}]
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
U10 = params{12};
U20 = params{13};
tol = params{14};
r = params{15};

% matrices of yn and zn values
yns = yn;
zns = zn;

% Initial energy, where potential energy is zero when the pendulum is at
% its lowest point.  Initially (yn,zn) = (y0,z0).
%E = (1/2)*norm(z0)^2 + g0*(l+y0(2)); 
E = (1/2)*norm(zn)^2 + g0*(l+yn(2)); 

for n=1:2
	% initial guesses for the simplified Newton method.
	if n==1
		x0F = [yn; zn; zn; U10];
		x0G = [zn; U20];
	else
		Yn1 = ynm;
		Yn2 = Y2;
		Zn1 = Z1;
		Zn2 = Z2;
		Un1 = U1;
		Un2 = U2;

		xnm = {ynm,znm};
		Xn = {Yn1,Yn2,Zn1,Zn2,Un1,Un2};
		SA_params = {ny,nz,nu,h,yn,zn,g0,l,b,r};

		Xnp0 = find_SA(xnm, Xn, SA_params);

		Ynp20 = Xnp0{1};
		Znp10 = Xnp0{2};
		Znp20 = Xnp0{3};
		znp0 = Xnp0{4};
		Unp10 = Xnp0{5};
		Unp20 = Xnp0{6};
		x0F = [Ynp20; Znp10; Znp20; Unp10];
		x0G = [znp0; Unp20];
	end

%%	%% solve the first nonlinear system
	% The function F needs to be redefined at every step using Fparams.
	% The argument of F is x = [Y2; Z1; Z2; U1].

	% Step size is h for the first step and rh for the second.
	% Note this changes Gparams as well.
	if n==1
		Fparams = {ny, nz, nu, h, yn, zn, g0, l, A, Ahat};
	else
		Fparams = {ny, nz, nu, r*h, yn, zn, g0, l, A, Ahat};
	end

	F = @(x) F1(x,Fparams);
	DF = @(x) DF1(x,Fparams);
	[F_root F_iter] = simplified_newton(F,DF(x0F),x0F,tol);
	%[F_root F_iter] = newton(F,DF,x0F,tol);

	% unpack F_root
	Y2 = F_root(1:ny);					% 1:2
	Z1 = F_root(ny+1:ny+nz);			% 3:4
	Z2 = F_root(ny+nz+1:ny+2*nz);		% 5:6
	U1 = F_root(ny+2*nz+1:ny+2*nz+nu);	% 7:7

%%	%% solve the second nonlinear system
	% The function G needs to be redefined at every step using Gparams.

	% Note that Gparams is Gparams with three new parameters.  This
	% makes the code easier to read and write but introduces some
	% inefficiency because not all of these parameters are actually
	% used.

	% The argument of G is x = [znp; U2].

	Gparams = [Fparams {b, Y2, U1}];
	G = @(x) G1(x,Gparams);
	DG = @(x) DG1(x,Gparams);
	[G_root G_iter] = simplified_newton(G,DG(x0G),x0G,tol);
	%[G_root G_iter] = newton(G,DG,x0G,tol);

	% unpack G_root 
	znp = G_root(1:nz);
	U2 = G_root(ny+1:ny+nu);

	% Update y
	ynp = Y2;

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
Y0_err = norm(Ynp20 - Y2);
Z0_err = norm([Znp10; Znp20] - [Z1; Z2]);
znp0_err = norm(znp0 - znp);
U0_err = norm([Unp10; Unp20] - [U1; U2]);

result = {Y0_err, Z0_err, znp0_err, U0_err, yns, zns, energy_vec};

