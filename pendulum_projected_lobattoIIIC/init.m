% init.m
% initial conditions for the pendulum example

g0 = 9.8; % Graviational constant.
alpha0 = [0;g0];
l = 1; % length of rod
h = 0.1; % stepsize
tol = 1e-6; % for the iterative solver (e.g. Newton's method)
N = 500; % number of steps

y0 = [3/5; -4/5]; % initial position
z0 = [-4/2; -3/2]; % initial velocity 

ny = 2; % dimensions of y
nz = 2; % dimensions of z
nu = 1; % dimensions of u

% Used in initial guesses for Newton's method
u0 = dot(z0,z0) - dot(y0,alpha0); % satisfies hidden constraint.  

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

% Note that params does NOT include h, N, y0, or z0.
params = {g0, l, tol, ny, nz, nu, method_str, b, c, A, Ahat, f, k, g, gyf};

% Uncomment the desired iterative scheme for the first nonlinear system
%iterative_scheme1 = 'newton';
%iterative_scheme1 = 'simplified newton 1';
iterative_scheme1 = 'simplified newton 2';

% Uncomment the desired iterative scheme for the second nonlinear system
%iterative_scheme2 = 'newton';
iterative_scheme2 = 'simplified newton';

