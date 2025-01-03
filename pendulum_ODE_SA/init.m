% init.m
% initial conditions and parameters for pendulum ODE example

%% Lobatto IIIA
%method_str = 'Lobatto IIIA';
%A = [0 0; 
%	1/2 1/2];
%b = [1/2; 1/2];
%c = [0; 1];

%% Lobatto IIIB
%method_str = 'Lobatto IIIB';
%A = [1/2 0; 
%	1/2 0];
%b = [1/2; 1/2];
%c = [0; 1];

%% Gauss
method_str = 'Gauss Method';
A = [1/4 1/4-sqrt(3)/6;
	1/4+sqrt(3)/6 1/4];
b = [1/2; 1/2];
c = [1/2-sqrt(3)/6; 1/2+sqrt(3)/6];

s = size(A,1); % number of stages

g0 = 9.8;
alpha0 = [0;g0];
l = 1; % length of the rod

y0 = [3/5; -4/5]; % initial position
z0 = [-4/2; -3/2]; % initial velocity
u0 = dot(z0,z0) - dot(y0,alpha0); % satisfies hidden constraint

psi0 = [y0; z0; u0]; % initial value
Psi0_trivial = kron(ones(s,1), psi0); % trivial starting approximation for first step

% dimensions of y, z, u
ny = size(y0,1);
nz = ny;
nu = 1;
npsi = ny+nz+nu;

% tolerance for Newton's method or simplified newton
tol = 1e-12; 


