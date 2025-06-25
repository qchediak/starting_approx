function dKdY_YZU = compute_dKdY(Y,Z,U,ny,nu,s)
	% function dKdY_YZU = compute_dKdY(Y,Z,U,ny,nu,s)
	% Returns the Jacobian matrix (dK/dY)(Y,Z,U), where Y=[Y_1;...;Y_s],
	% Z=[Z_1;...;Z_s], and U=[U1;...;Us], and
	% K(Y,Z,U)=[k(Y1,Z1,U1);...;k(Ys,Zs,Us)].
	%
	% This is specific to the pendulum problem.
	%
	% Inputs:
	%	Y: Column vector of internal stages [Y1;...;Ys]
	%	Z: Column vector of internal stages [Z1;...;Zs]
	%	U: Column vector of internal stages [U1;...;Us]
	%	ny: common dimension of Y1,...,Ys, 
	%	nu: common dimension of U1,...,Us
	% 	s: number of internal stages
	%
	% Outputs:
	%	dKdY_YZU: The (nu*ny*s)x(ny*s) matrix (dK/dY)(Y,Z,U).

	dKdY_YZU = (-1)*kron(blkdiag_vec(U,nu,s), eye(ny));
end

