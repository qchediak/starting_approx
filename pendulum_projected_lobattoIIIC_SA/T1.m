function v = T1(x,Tparams)
	% function v = T1(x,Tparams)
	% RHS of second linear system to solve for index 3 pendulum example
	% using projected Lobatto IIIC.
	%
	% Note that the linear system changes with each step, because (yn,zn)
	% are updated.  This is why we need to pass Tparams.
	%
	% Inputs:
	%	x = [znp; Unpsp]
	% 	Tparams = {ny, nz, nu, h, zn, b, k, gyf}
	%
	% Outputs:
	%	v: The value of T at x given Tparams.
	%
	% Note that Tparams has a cell array structure.

	% unpack Tparams
	ny = Tparams{1};
	nz = Tparams{2};
	nu = Tparams{3};
	h = Tparams{4};
	zn = Tparams{5};
	Ynp = Tparams{6};
	Znp = Tparams{7};
	Unp = Tparams{8};
	b = Tparams{9};
	k = Tparams{10};
	gyf = Tparams{11};

	s = size(b,1);

	% check that the dimensions of x are correct
	if size(x,1) ~= nz+nu
		error('Error (T1): x has the wrong number of rows.')
	elseif size(x,2) > 1
		error('Error (T1): x has too many columns.')
	end

	% unpack x
	znp = x(1:nz);
	Unpsp = x(nz+1:nz+nu);

	Ynpsp = Ynp((s-1)*ny+1:s*ny); % Y_{n+1,s+1} = Y_{n+1,s}
	Znpsp = Znp((s-1)*nz+1:s*nz); % Z_{n+1,s+1} = Z_{n+1,s}
	ynp = Ynpsp; % y_{n+1} = Y_{n+1,s}.  For readability.

	% Find T(x)
	v1 = znp - zn - h*kron([b(1:s-1)' 0],eye(nz)) * computeK(Ynp,Znp,Unp,k,ny,nu,s) - h*b(s)*k(Ynpsp,Znpsp,Unpsp);
	v2 = (1/h) * gyf(ynp,znp); % Note the scaling by 1/h

	v = [v1; v2];
end

