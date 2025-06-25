% init.m
% initial conditions for the pendulum example

g0 = 9.8; % Graviational constant.
alpha0 = [0;g0];
l = 1; % length of rod
tol = 1e-6; % for the iterative solver (e.g. Newton's method)

y0 = [3/5; -4/5]; % initial position
z0 = [-4/2; -3/2]; % initial velocity 

ny = 2; % dimensions of y
nz = 2; % dimensions of z
nu = 1; % dimensions of u

% Used in initial guesses for Newton's method
u0 = (1/l^2) * (dot(z0,z0) - dot(y0,alpha0)); % satisfies hidden constraint.  

% PRK coefficients: Uncomment the desired method
% Note that b is a column vector

%% Gauss
%method_str = 'Gauss';
%A = [1/4 1/4-sqrt(3)/6;
%	1/4+sqrt(3)/6 1/4];
%b = [1/2; 1/2];
%c = [1/2-sqrt(3)/6; 1/2+sqrt(3)/6];
%Ahat = A;

%% Radau IIA
method_str = 'Radau IIA';
A = [5/12 -1/12;
	3/4 1/4];
b = [3/4; 1/4];
c = [1/3; 1];
Ahat = A;

s = size(A,1); % number of stages

% RHS functions for the pendulum
f = @(y,z) z;
k = @(y,z,u) -u*y - alpha0;
g = @(y) dot(y,y) - l^2;

