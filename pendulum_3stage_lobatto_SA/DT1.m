function J = DT1(x,Tparams)
% function J = DT1(x,Tparams)
% Jacobian of the function T1.
%
% Note that the linear system changes with each step, because (yn,zn) are
% 	updated.  This is why we need to pass Tparams.
%
% Here 
% Tparams = {ny, nz, nu, h, yn, zn, g0, l, A0, A0hat, b, Ybar, Utilde}, 
% where h is stepsize; ny,nz,nu are respectively the dimensions of
% y,z,u; g0 is the gravitational constant; l is the length of the rod,
% A0 and A0hat are submatrices of the PRK coefficient matrices A and
% Ahat, b is the vector of weights of the PRK method, Ybar=[Y2; Y3],
% and Utilde = [U1; U2].
%
% Note that Tparams has a cell array structure.
%
% Here x = [znp; U3] are the variables being solved for.
%
% Note that this is specific to the pendulum problem.

% unpack Tparams
ny = Tparams{1};
nz = Tparams{2};
nu = Tparams{3};
h = Tparams{4};
yn = Tparams{5};
zn = Tparams{6};
g0 = Tparams{7};
l = Tparams{8};
A0 = Tparams{9};
A0hat = Tparams{10};
b = Tparams{11};
Ybar = Tparams{12};
Utilde = Tparams{13};

% check that the dimensions of x are correct
if size(x,1) ~= nz+nu
	error('Error (T1): x has the wrong number of rows.')
elseif size(x,2) > 1
	error('Error (T1): x has too many columns.')
end

% unpack x
znp = x(1:nz);
% U3 = x(nz+1:nz+nu); % not used

%Y1 = yn; % not used
%Y2 = Ybar(1:ny); % not used
Y3 = Ybar(ny+1:2*ny);
ynp = Y3;

% calculate the Jacobian as a 2x2 block matrix
J11 = eye(nz);
J12 = h*b(3)*Y3;
J21 = (1/h)*2*ynp';
J22 = 0;

J = [J11 J12; J21 J22];


