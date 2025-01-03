function J = computeDG(psi,params)
	% function DG_psi = computeDG(psi, params)
	% Jacobian of the RHS of the plane pendulum problem with the
	% constraint as an ODE.
	%
	% Inputs:
	%	psi: psi=(y,z,u), where y=(q_1,q_2), z=dy/dt
	%	params: cell array structure containing g0, ny, nz, nu
	%	g0: gravitational constant
	%	ny, nz, nu: dimension of y, z, and u, respectively 
	%		(should be 2, 2, 1 for this problem)
	%
	% Ouputs:
	%	J: This is DG(psi).

	% unpack params
	g0 = params{1};
	ny = params{2};
	nz = params{3};
	nu = params{4};

	% check dimensions of psi
	if size(psi,1) ~= ny+nz+nu
		error('Error (computeG): psi has the wrong number of rows.')
	elseif size(psi,2) > 1
		error('Error (computeG): psi is not a column vector')
	end

	% unpack psi
	y = psi(1:ny);
	z = psi(ny+1:ny+nz);
	u = psi(ny+nz+1:ny+nz+nu);

	alpha0 = [0; g0];

	% write DG(psi) as a 3x3 block matrix
	J11 = zeros(ny);
	J12 = eye(ny);
	J13 = zeros(ny,1);

	J21 = -u*eye(ny);
	J22 = zeros(ny);
	J23 = -y;

	J31 = (1 / dot(y,y)^2) * ( (8*u*dot(y,z) + 6*dot(z,alpha0)) * y' - (4*dot(y,y)*u) * z' );
	J32 = (-1 / dot(y,y)) * (4*u*y' + 3*alpha0');
	J33 = (-1 / dot(y,y)) * 4 * dot(y,z);

	J = [J11 J12 J13;
		J21 J22 J23;
		J31 J32 J33];

end


