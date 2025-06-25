function J = DG1(x,Gparams)
% function J = DG1(x,Gparams)
% Jacobian of the function G1.
%
% Note that the nonlinear system changes with each step, because (yn,zn) are
% 	updated.  This is why we need to pass Gparams.
%
% Here Gparams = {ny, nz, nu, h, yn, zn, g0, l, A, Ahat, b, Y2, U1}
% where h is stepsize; ny,nz,nu are respectively the dimensions of
% y,z,u; g0 is the gravitational constant; l is the length of the rod;
% A, Ahat, and b are the RK coefficient matrix; Y2 is the second
% internal stage of the y-component; and U1 is the first internal 
% stage of the u-component.
%
% Note that Gparams has a cell array structure.
%
% Here x = [znp; U2] are the variables being solved for.
%
% Note that this is specific to the pendulum problem.

% unpack Gparams
ny = Gparams{1};
nz = Gparams{2};
nu = Gparams{3};
h = Gparams{4};
yn = Gparams{5};
zn = Gparams{6};
g0 = Gparams{7};
l = Gparams{8};
A = Gparams{9};
Ahat = Gparams{10};
b = Gparams{11};
Y2 = Gparams{12};
U1 = Gparams{13};

% check that the dimensions of x are correct
if size(x,1) ~= nz+nu
	error('Error (G1): x has the wrong number of rows.')
elseif size(x,2) > 1
	error('Error (G1): x has too many columns.')
end

% unpack x
znp = x(1:nz);
U2 = x(nz+1:nz+nu);

alpha0 = [0; g0];
Y1 = yn;
ynp = Y2;

% calculate the Jacobian as a 2x2 block matrix
J11 = eye(2);
J12 = -h*b(2)*(-Y2);
J21 = 2*ynp';
J22 = 0;

J = [J11 J12; J21 J22];


