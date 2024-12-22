% 2024 Dec 10
% Solving the pendulum problem as an index 3 DAE using the 
% 3-stage lobatto method.
%
% With modifications this could be adapted to code which
% solves a general index 3 problem with n=m=2 and p=1 by a 
% PRK method with s=3.
%
% The purpose of this code is to solve the pendulum problem.  
%
% Note that this code does NOT use very good starting approximations.
%
% This assumes r=1.

clear; close all

init 

plots=true;
movie=true;

yn = y0;
zn = z0;

% matrices of yn and zn values
yns = yn;
zns = zn;

% Initial energy, where potential energy is zero when the pendulum is at
% its lowest point.
E = (1/2)*norm(z0)^2 + g0*(l+y0(2)); 

for n=1:N
	% initial guesses for the simplified Newton method.
	x0H = [yn; yn; zn; zn; zn; U10; U20];
	x0T = [zn; U20];

%%	%% solve the first nonlinear system
	% The function H needs to be redefined at every step using Hparams.
	% The argument of H is x = [Ybar; Z; Utilde], where
	% Ybar = [Y2; Y3], Z = [Z1; Z2; Z3], and Utilde = [U1; U2].


	Hparams = {ny, nz, nu, h, yn, zn, g0, l, A0, A0hat};
	H = @(x) H1(x,Hparams);
	DH = @(x) DH1(x,Hparams);
	[H_root H_iter] = simplified_newton(H,DH(x0H),x0H,tol);
	%[H_root H_iter] = newton(H,DH(x0H),x0H,tol);

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
	U1 = Utilde(1:nu);
	U2 = Utilde(nu+1:2*nu);

%%	%% solve the second nonlinear system
	% The function T needs to be redefined at every step using Tparams.

	% Note that Tparams is Tparams with three new parameters.  This
	% makes the code easier to read and write but not all of these
	% parameters are actually used.

	% The argument of T is x = [znp; U3].

	Tparams = [Hparams {b, Ybar, Utilde}];
	T = @(x) T1(x,Tparams);
	DT = @(x) DT1(x,Tparams);
	[T_root T_iter] = simplified_newton(T,DT(x0T),x0T,tol);
	%[T_root T_iter] = newton(T,DT,x0T,tol);

	% unpack T_root 
	znp = T_root(1:nz);
	U3 = T_root(nz+1:nz+nu);

	% Update y
	ynp = Y3;

	% update the matrix of yns, zns
	yns = [yns ynp];
	zns = [zns znp];

	% update energy
	Enp = (1/2)*norm(znp)^2 + g0*(l+ynp(2)); 
	E = [E Enp];

	zn = znp;
	yn = ynp;
end	

%% Plots
if plots 
	plots_params = {h,N};
	make_plots(yns,zns,E,plots_params)
end

if movie
	make_movie(yns,zns,plots_params)
end

