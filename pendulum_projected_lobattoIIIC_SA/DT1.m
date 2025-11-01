function J = DT1(x,Tparams)
	% function J = DT1(x,Tparams)
	% Jacobian of the function T1.
	%
	% Note that the linear system changes with each step, because (yn,zn) are
	% 	updated.  This is why we need to pass Tparams.
	%
	% Inputs:
	%	x = [znp; Unpsp]
	% 	Tparams = {ny, nz, nu, h, zn, b, k, gyf}
	%
	% Outputs:
	%	J: The value of DT at x given Tparams.

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

	ynp = Ynp((s-1)*ny+1:s*ny); % y_{n+1} = Y_{n+1,s}
	

	% calculate the Jacobian as a 2x2 block matrix
	J11 = eye(nz);
	J12 = h*b(s)*ynp;
	J21 = (1/h)*2*ynp'; % scale by 1/h so the Jacobian remains invertible
	J22 = 0;

	J = [J11 J12; J21 J22];
end

