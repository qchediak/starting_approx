function G_psi = computeG(psi, params)
	% function G_psi = computeG(psi, params)
	% RHS of the plane pendulum problem with the constraint as an ODE
	%
	% Inputs:
	%	psi: psi=(y,z,u), where y=(q_1,q_2), z=dy/dt
	%	params: cell array structure containing g0, ny, nz, nu
	%	g0: gravitational constant
	%	ny, nz, nu: dimension of y, z, and u, respectively 
	%		(should be 2, 2, 1 for this problem)
	%
	% Ouputs:
	%	G_psi: This is G(psi).

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

	v1 = z;
	v2 = -u*y - alpha0;
	v3 = ( 1 / dot(y,y) ) * (-4*u*dot(y,z) -3*dot(z,alpha0));

	G_psi = [v1; v2; v3];

end

