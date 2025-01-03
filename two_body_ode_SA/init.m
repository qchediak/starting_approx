% init.m
% initial conditions and parameters for two-body ODE example

%% Parameters relating to the numerical method

% uncomment the desired method
%method_str = 'Lobatto IIIA';
%method_str = 'Lobatto IIIB';
method_str = 'Gauss Method';

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
tol = 1e-12; 

%% Parameters relating to the physical problem
d = 2; % common dimension of q and p
e = 0.6; % Eccentricity of the orbit.  Only used in initial values.

%% starting values 
% from Geometric Numerical Integration by Hairer et al, page 12
q0 = [1-e; 0];
p0 = [0; sqrt((1+e)/(1-e))];

y0 = [q0; p0];

Y0_trivial = kron(ones(s,1), y0); % trivial starting approximation

