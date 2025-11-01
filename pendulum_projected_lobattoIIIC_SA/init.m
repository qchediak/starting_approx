% init.m
% initial conditions for the pendulum example

g0 = 9.8; % Graviational constant.
alpha0 = [0;g0];
l = 1; % length of rod
tol = 1e-12; % for the iterative solver (e.g. Newton's method)



y0 = [3/5; -4/5]; % initial position
z0 = [-4/2; -3/2]; % initial velocity 

ny = 2; % dimensions of y
nz = 2; % dimensions of z
nu = 1; % dimensions of u

% Used in initial guesses for Newton's method
u0 = (1/l^2) * (dot(z0,z0) - dot(y0,alpha0)); % satisfies hidden constraint.  

% PRK coefficients: Uncomment the desired method

% 3-stage Lobatto IIIC PRK coefficients
method_str = 'Lobatto IIIC';
A = [1/6 -1/3 1/6;
	1/6 5/12 -1/12;
	1/6 2/3 1/6];
b = [1/6; 2/3; 1/6];
c = [0; 1/2; 1];
Ahat = A;

s = size(A,1); % number of stages

% RHS functions for the pendulum
f = @(y,z) z;
k = @(y,z,u) -u*y - alpha0;
g = @(y) dot(y,y) - l^2;
gyf = @(y,z) 2*dot(y,z);

% Uncomment the desired iterative scheme for the first nonlinear system
%iterative_scheme1 = 'newton';
%iterative_scheme1 = 'simplified newton 1';
iterative_scheme1 = 'simplified newton 2';

% Uncomment the desired iterative scheme for the second nonlinear system
%iterative_scheme2 = 'newton';
iterative_scheme2 = 'simplified newton';

% Options for plotting errors:
% 1: Useful for testing.  Allows more flexibility in what plots to make,
% 	such as component or composite (see below), and separately plots
% 	znp0 and Unpsp0 errors.
% 2: The option used to make the graphs in my thesis.  This option plots
% 	the y-components separately, combines znp0 with Z0, and Unpsp0 with
% 	U0.  This option assumes the method is projected Lobatto IIIC and 
%	s=3.
error_plotting_option = 2;

% Options for plotting errors: 
% Either 'composite' (i.e. error in all of X_{n+1}^{(0)})
% or 'component' (i.e. error in each X_{n+1,i}^{(0)} for i=1,...,s).
% Currently the 'component' option is only available for s=3.
error_plotting_Y = 'component';
error_plotting_Z = 'component';
error_plotting_U = 'component';

