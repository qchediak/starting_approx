function J = DF1(x,Fparams)
% function J = DF1(x,Fparams)
% Jacobian of the function F1.
%
% Here x = [Y2; Z1; Z2; U1] are the variables being solved for and 
% Fparams = {ny, nz, nu, h, yn, zn, g0, l, A, Ahat} are the parameters needed.
% Note that Fparams has a cell array structure.
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
	error('Error (DF1): x has the wrong number of rows.')
elseif size(x,2) > 1
	error('Error (DF1): x has too many columns.')
end

% unpack x
Y2 = x(1:ny);
Z1 = x(ny+1:ny+nz);	
Z2 = x(ny+nz+1:ny+2*nz);
U1 = x(ny+2*nz+1:ny+2*nz+nu);

Y1 = yn;
alpha0 = [0; g0];

% calculate the Jacobian as a 4x4 block matrix
J11 = eye(2);
J12 = -h*A(2,1)*eye(2);
J13 = -h*A(2,2)*eye(2);
J14 = zeros(2,1);

J21 = zeros(2,2);
J22 = eye(2);
J23 = zeros(2,2);
J24 = -h*Ahat(1,1)*(-Y1);

J31 = zeros(2,2);
J32 = zeros(2,2);
J33 = eye(2);
J34 = -h*Ahat(2,1)*(-Y1);

J41 = (1/h^2)*2*Y2'; % multiply by 1/h^2 so Jacobian is invertible
J42 = zeros(1,2);
J43 = zeros(1,2);
J44 = 0;

J = [J11 J12 J13 J14; 
	J21 J22 J23 J24; 
	J31 J32 J33 J34;
	J41 J42 J43 J44];

