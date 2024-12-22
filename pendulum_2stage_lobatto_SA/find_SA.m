function Xnp0 = find_SA(xnm, Xn, SA_params)
% function Xnp0 = find_SA(xnm, Xn, SA_params)
% Finds the starting approximations Xnp0 to Xnp based on the previous 
% internal stages Xn and the past step xnm.
%
% Inputs:
%	xnm = {ynm, znm} 
% 	Xn = {Yn1,Yn2,Zn1,Zn2,Un1,Un2}
%	SA_params = {ny,nz,nu,h,yn,zn,g0,l,b,r}
%	ny, nz, nu: dimensions of y,z,u
%	h: stepsize of the first step (does NOT factor in r)
%	(yn, zn): The first step after (ynm,znm)
%	g0: gravitational constant near Earth's surface (9.8 m/s^2)
%	l: length of rod
%	b: weights of the PRK method
%	r: ratio of second stepsize / first stepsize
%
% Outputs:
% 	Xnp0 = {Ynp20, Znp10, Znp20, znp0, Unp10, Unp20}
%
% Note that xnm, Xn, params, and Xnp0 all have a cell array structure.
%
% Assumes we are using Lobatto IIIA-IIIB with s=2 to solve the plane
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
b = SA_params{9};
r = SA_params{10};

% unpack xnm
ynm = xnm{1};
znm = xnm{2};

% unpack Xn
Yn1 = Xn{1};
Yn2 = Xn{2};
Zn1 = Xn{3};
Zn2 = Xn{4};
Un1 = Xn{5};
Un2 = Xn{6};

% RHS functions.  Note that these depend on g0 and l.
f = @(y,z,u) z;
k = @(y,z,u) -u*y - [0; g0];
%g = @(y) y' * y - l^2; % not used

beta11_star = 0; % only needed for matrix calculation
beta12_star = 0; % only needed for matrix calculation
beta21_star = -r;
beta22_star = 0;

beta11_hat_star = -(1/2) * r;
beta12_hat_star = 0;
beta21_hat_star = -(1/2) * r;
beta22_hat_star = 0;
beta31_hat_star = -r;
beta32_hat_star = 0;

D11_star = 1/2;
D12_star = 1/2;
D21_star = 1/2;
D22_star = 1/2;

beta_star = [beta11_star beta12_star; 
				beta21_star beta22_star];

beta_hat_star = [beta11_hat_star beta12_hat_star; 
					beta21_hat_star beta22_hat_star;
					beta31_hat_star beta32_hat_star];

D_star = [D11_star D12_star; 
			D21_star D22_star];

%% Find the unstarred coefficients
E = [0 1; 1 0];
ones2 = ones(2,1);
ones3 = ones(3,1);

beta = kron(ones2,b') - beta_star * E;
%beta_hat = kron(ones2,b') - beta_hat_star * E;
beta_hat = kron(ones3,b') - beta_hat_star * E;

%beta31_hat = b(1) - beta32_hat_star;
%beta32_hat = b(2) - beta31_hat_star;

D = D_star * E;

%% Find the starting approximations Xnp0

F_n = [f(Yn1,Zn1); f(Yn2,Zn2)];
K_n = [k(Yn1,Zn1,Un1); k(Yn2,Zn2,Un2)];

Ynp0 = kron(ones2,ynm) + h * kron(beta,eye(ny)) * F_n;
%Znp0 = kron(ones2,znm) + h * kron(beta_hat,eye(nz)) * K_n;
%znp0 = znm + h * ( beta31_hat * k(Yn1,Zn1,Un1) + beta32_hat * k(Yn2,Zn2,Un2) );
Znp0 = kron(ones3,znm) + h * kron(beta_hat,eye(nz)) * K_n;
Unp0 = kron(D,eye(nu)) * [Un1; Un2];

% unpack Ynp0, Znp0, Unp0
Ynp10 = Ynp0(1:ny); % not needed
Ynp20 = Ynp0(ny+1:2*ny);
Znp10 = Znp0(1:nz);
Znp20 = Znp0(nz+1:2*nz);
znp0 = Znp0(2*nz+1:3*nz); % Znp30
Unp10 = Unp0(1:nu);
Unp20 = Unp0(nu+1:2*nu);

Xnp0 = {Ynp20,Znp10,Znp20,znp0,Unp10,Unp20};


