% 2024 Nov 27
% Solving the pendulum problem as an index 3 DAE using the 
% 2-stage lobatto method.
%
% With minor modifications this could be adapted to code which
% solves a general index 3 problem with n=m=2 and p=1 by a 
% PRK method with s=2.
%
% The purpose of this code is to solve the pendulum problem.  
%
% Note that this code does NOT use very good starting approximations.

clear; close all

init 

plots=true;
movie=false;

yn = y0;
zn = z0;

% matrices of yn and zn values
yns = yn;
zns = zn;

% initial energy, where potential energy is zero when the pendulum is at
% its lowest point.
E = (1/2)*norm(z0)^2 + g0*(l+y0(2)); 

for n=1:N
	% initial guesses for the simplified Newton method.
	x0F = [yn; zn; zn; U10];
	x0G = [zn; U20];

%%	%% solve the first nonlinear system
	% The function F needs to be redefined at every step using Fparams.
	% The argument of F is x = [Y2; Z1; Z2; U1].

	Fparams = {ny, nz, nu, h, yn, zn, g0, l, A, Ahat};
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

	zn = znp;
	yn = ynp;
end	

%% Plots
if plots 
	plots_params = {h,N,s};
	make_plots(yns,zns,E,plots_params)
end

if movie
	make_movie(yns,zns,plots_params)
end


