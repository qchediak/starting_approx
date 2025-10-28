function v = F1(x,Fparams)
% function v = F1(x,Fparams)
%
% LHS of first nonlinear system to solve for index 3 pendulum example.
%
% Note that the nonlinear system changes with each step, because (yn,zn) are
% 	updated.  This is why we need to pass Fparams.
%
% Here Fparams = {ny, nz, nu, h, yn, zn, g0, l, A, Ahat}, 
% where h is stepsize; ny,nz,nu are respectively the dimensions of
% y,z,u; g0 is the gravitational constant; l is the length of the rod,
% and A and Ahat are the PRK coefficient matrices.  
%
% Note that Fparams has a cell array structure.
%
% Here x = [Y2; Z1; Z2; U1] are the variables being
% solved for.  
%
% Note that this is specific to the pendulum problem.

% unpack Fparams
ny = Fparams{1};
nz = Fparams{2};
nu = Fparams{3};
h = Fparams{4};
yn = Fparams{5};
zn = Fparams{6};
g0 = Fparams{7};
l = Fparams{8};
A = Fparams{9};
Ahat = Fparams{10};

% check that the dimensions of x are correct
if size(x,1) ~= ny+2*nz+nu
	error('Error (F1): x has the wrong number of rows.')
elseif size(x,2) > 1
	error('Error (F1): x has too many columns.')
end

% unpack x
Y2 = x(1:ny);
Z1 = x(ny+1:ny+nz);
Z2 = x(ny+nz+1:ny+2*nz);
U1 = x(ny+2*nz+1:ny+2*nz+nu);

alpha0 = [0; g0];
Y1 = yn;

% calculate F(x) with above parameters
v1 = Y2 - yn - h*(A(2,1)*Z1 + A(2,2)*Z2);
v2 = Z1 - zn - h*Ahat(1,1)*(-U1*Y1 - alpha0);
v3 = Z2 - zn - h*Ahat(2,1)*(-U1*Y1 - alpha0);
v4 = (1/h^2) * (dot(Y2,Y2) - l^2);

v = [v1; v2; v3; v4];

