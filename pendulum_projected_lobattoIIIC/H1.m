function v = H1(x,Hparams)
% function v = H1(x,Hparams)
%
% The nonlinear system to solve for index 3 pendulum example.
%
% Note that the linear system changes with each step, because (yn,zn) are
% 	updated.  This is why we need to pass Hparams.
%
% Here Hparams = {ny, nz, nu, h, yn, zn, g0, l, A, Ahat}, 
% where h is stepsize; ny,nz,nu are respectively the dimensions of
% y,z,u; g0 is the gravitational constant; l is the length of the rod,
% and A and Ahat are the PRK coefficient matrices.  
%
% Note that Hparams has a cell array structure.
%
% Here x = [Y; Z; U] are the variables being
% solved for.  
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
alpha0 = [0;g0];

% check that the dimensions of x are correct
if size(x,1) ~= s*(ny+nz+nu)
	error('Error (H1): x has the wrong number of rows.')
elseif size(x,2) > 1
	error('Error (H1): x has too many columns.')
end

% unpack x
Y = x(1:ny*s);				% for s=2, ny=nz=2, nu=1, this is 1:4
Z = x(ny*s+1:ny*s+nz*s);	% for s=2, ny=nz=2, nu=1, this is 5:8
U = x(ny*s+nz*s+1:end);		% for s=2, ny=nz=2, nu=1, this is 9:10

alpha0 = [0; g0];

ones_vec = ones(s,1);

% calculate H(x) with above parameters
v1 = Y - kron(ones_vec, yn) - h * kron(A, eye(ny)) * computeF(Y,Z,f,ny,s);
v2 = Z - kron(ones_vec, zn) - h * kron(Ahat, eye(nz)) * computeK(Y,Z,U,k,ny,nu,s);
v3 = computeG(Y,g,ny,s);
%v3 = (1/h^2) * computeG(Y,g,ny,s); % scaled constraint

v = [v1; v2; v3];

