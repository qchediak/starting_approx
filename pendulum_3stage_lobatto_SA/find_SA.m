function Xnp0 = find_SA(xnm, Xn, SA_params)
% function Xnp0 = find_SA(xnm, Xn, SA_params)
% Finds the starting approximations Xnp0 to Xnp based on the previous 
% internal stages Xn and the past step xnm.
%
% Inputs:
%	xnm = {ynm, znm} 
% 	Xn = {Yn1,Yn2,Yn3,Zn1,Zn2,Zn3,Un1,Un2,Un3}
%	SA_params = {ny,nz,nu,h,yn,zn,g0,l,A,Ahat,b,c,r}
%	ny, nz, nu: dimensions of y,z,u
%	h: stepsize of the first step (does NOT factor in r)
%	(yn, zn): The first step after (ynm,znm)
%	g0: gravitational constant near Earth's surface (9.8 m/s^2)
%	l: length of rod
%	A, Ahat: PRK coefficient matrices
%	b: weights of the PRK method
%	r: ratio of second stepsize / first stepsize
%
% Outputs:
% 	Xnp0 = {Ynp20, Ynp30, Znp10, Znp20, Znp30, znp0, Unp10, Unp20, Unp30}
%
% Note that xnm, Xn, params, and Xnp0 all have a cell array structure.
%
% Assumes we are using Lobatto IIIA-IIIB with s=3 to solve the plane
% pendulum formulated as an index 3 DAE.

% unpack params
ny = SA_params{1};
nz = SA_params{2};
nu = SA_params{3};
h = SA_params{4};
yn = SA_params{5};
zn = SA_params{6};
g0 = SA_params{7};
l = SA_params{8};
A = SA_params{9};
Ahat = SA_params{10};
b = SA_params{11};
c = SA_params{12};
r = SA_params{13};

% unpack xnm
ynm = xnm{1};
znm = xnm{2};

% unpack Xn
Yn1 = Xn{1};
Yn2 = Xn{2};
Yn3 = Xn{3};
Zn1 = Xn{4};
Zn2 = Xn{5};
Zn3 = Xn{6};
Un1 = Xn{7};
Un2 = Xn{8};
Un3 = Xn{9};

% RHS functions.  Note that these depend on g0 and l.
f = @(y,z,u) z;
k = @(y,z,u) -u*y - [0; g0];
%g = @(y) y' * y - l^2; % not used

%% find starred coefficients
beta_star = [0 0 0;
			-r*c(2)*(1+r*c(2)) r^2*c(2)^2 0;
			-r*c(3)*(1+r*c(3)) r^2*c(3)^2 0];

beta_hat_star = [-r*c(1)-r^2*Ahat(1,2) r^2*Ahat(1,2) 0;
					-r*c(2)-r^2*Ahat(2,2) r^2*Ahat(2,2) 0;
					-r*c(3)-r^2*Ahat(3,2) r^2*Ahat(3,2) 0;
					-r^2-r r^2 0];

D_star = [1+2*r*c(1) -2*r*c(1) 0;
			1+2*r*c(2) -2*r*c(2) 0;
			1+2*r*c(3) -2*r*c(3) 0];

%% find unstarred coefficients

E = [0 0 1; 
	0 1 0; 
	1 0 0];

beta = kron(ones(3,1),b') - beta_star*E;
beta_hat = kron(ones(4,1),b') - beta_hat_star*E;
D = D_star*E;

F = [f(Yn1,Zn1); f(Yn2,Zn2); f(Yn3,Zn3)];
K = [k(Yn1,Zn1,Un1); k(Yn2,Zn2,Un2); k(Yn3,Zn3,Un3)];
Un = [Un1; Un2; Un3];

%% find the starting approximations using tensor notation

Ynp0 = kron(ones(3,1),ynm) + h * kron(beta, eye(ny))*F;
Znp0 = kron(ones(4,1),znm) + h * kron(beta_hat, eye(nz))*K;
Unp0 = kron(D,eye(nu))*Un;

% unpack Ynp0, Znp0, Unp0

% Ynp10 = Ynp0(1:ny); % not used
Ynp20 = Ynp0(ny+1:2*ny);
Ynp30 = Ynp0(2*ny+1:3*ny);
Znp10 = Znp0(1:nz);
Znp20 = Znp0(nz+1:2*nz);
Znp30 = Znp0(2*nz+1:3*nz);
znp0 = Znp0(3*nz+1:4*nz); % Znp40
Unp10 = Unp0(1:nu);
Unp20 = Unp0(nu+1:2*nu);
Unp30 = Unp0(2*nu+1:3*nu);

Xnp0 = {Ynp20,Ynp30,Znp10,Znp20,Znp30,znp0,Unp10,Unp20,Unp30};


