function J = DH1(x,Hparams)
% function J = DH1(x,Hparams)
% Jacobian of the function H1.
%
% Note that the linear system changes with each step, because (yn,zn) are
% updated.  This is why we need to pass Hparams.
%
% Here Hparams = {ny, nz, nu, h, yn, zn, g0, l, A0, A0hat}, 
% where h is stepsize; ny,nz,nu are respectively the dimensions of
% y,z,u; g0 is the gravitational constant; l is the length of the rod,
% and A0 and A0hat are submatrices of the PRK coefficient matrices A
% and Ahat.
%
% Note that Hparams has a cell array structure.
%
% Here x = [Ybar; Z; Utilde] are the variables being solved for,
% where Ybar = [Y2; Y3], Z=[Z1; Z2; Z3], and Utilde=[U1; U2].
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
A0 = Hparams{9};
A0hat = Hparams{10};

% check that the dimensions of x are correct
if size(x,1) ~= 2*ny+3*nz+2*nu
	error('Error (H1): x has the wrong number of rows.')
elseif size(x,2) > 1
	error('Error (H1): x has too many columns.')
end

% unpack x
Ybar = x(1:2*ny);
Z = x(2*ny+1:2*ny+3*nz);
Utilde = x(2*ny+3*nz+1:2*ny+3*nz+2*nu);

% unpack Ybar
Y2 = Ybar(1:ny);
Y3 = Ybar(ny+1:2*ny);

% unpack Utilde
U1 = Utilde(1:nu);
U2 = Utilde(nu+1:2*nu);

alpha0 = [0; g0];
Y1 = yn;
ones2 = ones(2,1);
ones3 = ones(3,1);

% Define F, Ktilde, G, 
% where F = [f(Y1,Z1); f(Y2,Z2); f(Y3,Z3)], 
% Ktilde = [k(Y1,Z1,U1); k(Y2,Z2,U2)], 
% and G = [g(Y2); g(Y3)].
F = Z;
Ktilde = [-U1*Y1-alpha0; -U2*Y2-alpha0];
G = [dot(Y2,Y2)-l^2; dot(Y3,Y3)-l^2];

% calculate the Jacobian as a 3x3 block matrix
J11 = eye(2*ny);
J12 = -h*kron(A0,eye(ny));
J13 = zeros(2*ny,2);

J21 = -h*kron(A0hat,eye(nz))*[zeros(ny) zeros(ny); -U2*eye(ny) zeros(ny)];
J22 = eye(3*nz);
J23 = -h*kron(A0hat,eye(nz))*[-Y1 zeros(2,1); zeros(2,1) -Y2];

J31 = (1/h^2)*[2*Y2' zeros(1,2); zeros(1,2) 2*Y3']; % scale by 1/h^2 so the Jacobian remains invertible
J32 = zeros(2,3*nz);
J33 = zeros(2,2);

J = [J11 J12 J13; 
	J21 J22 J23; 
	J31 J32 J33];

