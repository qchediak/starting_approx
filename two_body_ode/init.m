% init.m
% initial conditions and parameters for two-body ODE example

%% Parameters relating to the physical problem
d = 2; % common dimension of q and p
e = 0.6; % eccentricity of the orbit

%% starting values 
% from Geometric Numerical Integration by Hairer et al, page 12
q0 = [1-e; 0];
p0 = [0; sqrt((1+e)/(1-e))];

y0 = [q0; p0];

E0 = (1/2)*norm(p0)^2 - 1/norm(q0); % initial energy

% Semi-major axis of elliptical orbit. 
% See, for example, "Classical Dynamics of Particles and Systems",
% Fifth Edition, by Thornton and Marion, eq. 8.42, p. 301.
% (The constant k=GMm in that equation is equal to 1 for our problem.)
a = 1 / (2*abs(E0)); 

% A constant used in the plot of the true solution trajectory.
% Ibid, eq. 8.42, p. 301
alpha0 = a*(1-e^2);

% Orbital period in terms of a.
% See, for example, "Classical Dynamics of Particles and Systems",
% Fifth Edition, by Thornton and Marion, eq. 8.49, p. 303.
% (The expression GM in that equation is equal to 1 for our problem.)
T = 2*pi*a^(3/2); 
 
%% Parameters relating to the numerical method

h = 0.1;
N = round(8*T/h); % choose N so that 8 orbital periods are shown

% uncomment the desired method
method_str = 'Lobatto IIIA';
%method_str = 'Lobatto IIIB';
%method_str = 'Gauss Method';

if strcmp(lower(strrep(method_str, ' ', '')), 'lobattoiiia')
	% Lobatto IIIA
	A = [0 0; 
		1/2 1/2];
	b = [1/2; 1/2];
	c = [0; 1];
elseif strcmp(lower(strrep(method_str, ' ', '')), 'lobattoiiib')
	% Lobatto IIIB
	A = [1/2 0; 
		1/2 0];
	b = [1/2; 1/2];
	c = [0; 1];
elseif strcmp(lower(strrep(method_str, ' ', '')), 'gaussmethod') || strcmp(lower(strrep(method_str, ' ', '')), 'gauss')
	% Gauss
	method_str = 'Gauss Method';
	A = [1/4 1/4-sqrt(3)/6;
		1/4+sqrt(3)/6 1/4];
	b = [1/2; 1/2];
	c = [1/2-sqrt(3)/6; 1/2+sqrt(3)/6];
else
	fprintf('Error (init): Invalid value of method_str')
end

s = size(A,1); % number of stages

% tolerance for Newton's method or simplified newton
tol = 1e-10; 

