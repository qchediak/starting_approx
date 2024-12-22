% init.m
% initial conditions for the pendulum example

g0 = 9.8; % Graviational constant.
l = 1; % length of rod
tol = 1e-6; % for Newton's method

y0 = [3/5; -4/5]; % initial position
z0 = [-4/2; -3/2]; % initial velocity 

ny = 2; % dimensions of y
nz = 2; % dimensions of z
nu = 1; % dimensions of u

% Used as the initial guess for Newton's method in the first step
U10 = 10;
U20 = 10;

% RK coefficients of Lobatto IIIA-IIIB with s=2
% Note that b is a column vector
a11 = 0;
a12 = 0;
a21 = 1/2;
a22 = 1/2;

a11hat = 1/2;
a12hat = 0;
a21hat = 1/2;
a22hat = 0;

b1 = 1/2;
b2 = 1/2;

c1 = 0;
c2 = 1;

A = [a11 a12; a21 a22];
Ahat = [a11hat a12hat; a21hat a22hat];
b = [b1; b2]; 
c = [c1; c2];

yn=y0;
zn=z0;


