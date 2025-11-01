function J = computeJ(un,Hparams)
	% function J = computeJ(un,Hparams)
	% Finds the matrix J for use in the root finder iterative_dae.
	%
	% Specific to the pendulum problem.
	%
	% Inputs:
	%	un: The value u_n such that (y_n, z_n, u_n) satisfies the hidden
	%		constraint
	%	Hparams: listed below
	%
	% Outputs:
	%	J: The matrix J used in the iterative root finder.

	% unpack Hparams
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

	% make J as a 3x3 block matrix

	J11 = eye(ny*s);
	J12 = -h * kron(A,eye(ny));
	J13 = zeros(ny*s,nu*s);

	%J21 = -h * kron(Ahat, eye(nz)) * compute_dKdY(Y,Z,U,ny,nu,s);
	J21 = -h * kron(Ahat, eye(nz)) * compute_dKdY(kron(ones(s,1),yn),kron(ones(s,1),zn),kron(ones(s,1),un),ny,nu,s);
	J22 = eye(ny*s);
	%J23 = -h * kron(Ahat, eye(nz)) * compute_dKdU(Y,Z,U,ny,s);
	J23 = -h * kron(Ahat, eye(nz)) * compute_dKdU(kron(ones(s,1),yn),kron(ones(s,1),zn),kron(ones(s,1),un),ny,s);

	J31 = zeros(nu*s,ny*s);
	J32 = compute_dGdY(kron(ones(s,1),yn),ny,s) * kron(A,eye(ny));
	J33 = zeros(nu*s,nu*s);

	J = [J11 J12 J13;
		J21 J22 J23;
		J31 J32 J33];


