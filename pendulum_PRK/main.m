% 2025 June 24
% Solving the pendulum problem as an index 3 DAE using a 
% PRK method.
%
% The purpose of this code is to solve the pendulum problem.  
% It is not intended as a general index 3 DAE solver.
%
% Note that this code uses only trivial starting approximations for y
% and z.

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
	% This is the value of un which satisfies the hidden constraint.
	% Used in the initial guess for the iterative solver.
	un = (1/l^2) * (dot(zn,zn) - dot(yn,alpha0)); 

	% initial guesses for the simplified Newton method.
	x0 = [kron(ones(s,1),yn); kron(ones(s,1),zn); kron(ones(s,1),un)];

	[ynp znp Ynp Znp Unp iter] = solve_PRK_DAE_single_step(params, h, yn, zn, un, x0, 'simplified newton 2');

	% update the matrix of yns, zns
	yns = [yns ynp];
	zns = [zns znp];

	% update energy
	Enp = (1/2)*norm(znp)^2 + g0*(l+ynp(2)); 
	E = [E Enp];

	zn = znp;
	yn = ynp;
end	

plots_params = {h,N,s,method_str};

%% Plots
if plots 
	make_plots(yns,zns,E,plots_params)
end

if movie
	make_movie(yns,zns,plots_params)
end


