function J = DH1(x,Hparams)
% function J = DH1(x,Hparams)
% Jacobian of the function H1.
%
% Here x = [Y; Z; U] are the variables being solved for and 
% Hparams = {ny, nz, nu, h, yn, zn, g0, l, A, Ahat} are the parameters needed.
% Note that Hparams has a cell array structure.
%
% Note that this is specific to the pendulum problem.

% unpack Hparams
ny = Hparams{1};
nz = Hparams{2};
nu = Hparams{3};
h = Hparams{4};
yn = Hparams{5};
zn = Hparams{6};
g0 = Hparams{7};
l = Hparams{8};
A = Hparams{9};
Ahat = Hparams{10};
f = Hparams{11};
k = Hparams{12};
g = Hparams{13};

s = size(A,1);

% check that the dimensions of x are correct
if size(x,1) ~= s*(ny+nz+nu)
	error('Error (DH1): x has the wrong number of rows.')
elseif size(x,2) > 1
	error('Error (DH1): x has too many columns.')
end

% unpack x
Y = x(1:ny*s);				% for s=2, ny=nz=2, nu=1, this is 1:4
Z = x(ny*s+1:ny*s+nz*s);	% for s=2, ny=nz=2, nu=1, this is 5:8
U = x(ny*s+nz*s+1:end);		% for s=2, ny=nz=2, nu=1, this is 9:10

% compute the derivatives of K at (Y,Z,U)
dKdY = compute_dKdY(Y,Z,U,ny,nu,s);
dKdU = compute_dKdU(Y,Z,U,ny,s);
dGdY = compute_dGdY(Y,ny,s);
%dGdY = (1/h^2) * compute_dGdY(Y,ny,s); % scaled constraint

% calculate the Jacobian as a 3x3 block matrix
J11 = eye(ny*s);
J12 = -h * kron(A, eye(ny));
J13 = zeros(ny*s,s);

J21 = -h * kron(Ahat, eye(ny)) * dKdY;
J22 = eye(ny*s);
J23 = -h * kron(Ahat, eye(ny)) * dKdU;

%J31 = compute_dGdY(Y,ny,s);
J31 = dGdY;
J32 = zeros(s,ny*s);
J33 = zeros(s,s);

J = [J11 J12 J13; 
	J21 J22 J23; 
	J31 J32 J33];


