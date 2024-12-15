% init.m
% initial conditions for the pendulum example

g0 = 9.8; % Graviational constant.
l = 1; % length of rod
h = 0.1; % stepsize
tol = 1e-6; % for Newton's method
N = 500; % number of steps

y0 = [3/5; -4/5]; % initial position
z0 = [-4/2; -3/2]; % initial velocity 

ny = 2; % dimensions of y
nz = 2; % dimensions of z
nu = 1; % dimensions of u

% Used in initial guesses for Newton's method
U10 = 0;
U20 = 0;
U30 = 0;

% PRK coefficients of Lobatto IIIA-IIIB with s=3

A = [0 0 0;
	5/24 1/3 -1/24;
	1/6 2/3 1/6];

Ahat = [1/6 -1/6 0;
		1/6 1/3 0;
		1/6 5/6 0];

b = [1/6; 2/3; 1/6];
c = [0; 1/2; 1];

A0 = A(2:3,:);
A0hat = Ahat(:,1:2);


