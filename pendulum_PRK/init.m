% init.m
% initial conditions for the pendulum example

g0 = 9.8; % Graviational constant.
alpha0 = [0;g0];
l = 1; % length of rod
h = 5e-3; % stepsize
tol = 1e-6; % for the iterative solver (e.g. Newton's method)
N = 500; % number of steps

y0 = [3/5; -4/5]; % initial position
z0 = [-4/2; -3/2]; % initial velocity 

ny = 2; % dimensions of y
nz = 2; % dimensions of z
nu = 1; % dimensions of u

% Used in initial guesses for Newton's method
u0 = dot(z0,z0) - dot(y0,alpha0); % satisfies hidden constraint.  

% PRK coefficients: Uncomment the desired method
% Note that b is a column vector

%% 3-stage Gauss
method_str = 'Gauss Method';
A=[5/36, 2/9-sqrt(15)/15, 5/36-sqrt(15)/30;
	5/36+sqrt(15)/24, 2/9, 5/36-sqrt(15)/24;
	5/36+sqrt(15)/30, 2/9+sqrt(15)/15, 5/36];
b = [5/18; 4/9; 5/18];
c = [1/2-sqrt(15)/10; 1/2; 1/2+sqrt(15)/10];
Ahat = A;

%% 2-stage Radau IIA
%method_str = 'Radau IIA';
%A = [5/12 -1/12;
%	3/4 1/4];
%b = [3/4; 1/4];
%c = [1/3; 1];
%Ahat = A;

s = size(A,1); % number of stages

% RHS functions for the pendulum
f = @(y,z) z;
k = @(y,z,u) -u*y - alpha0;
g = @(y) dot(y,y) - l^2;

% Note that params does NOT include h, N, y0, or z0.
params = {g0, l, tol, ny, nz, nu, method_str, b, c, A, Ahat, f, k, g};

