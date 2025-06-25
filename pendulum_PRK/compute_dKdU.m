function dKdU_YZU = compute_dKdU(Y,Z,U,ny,s)
	% function dKdU_YZU = compute_dKdU(Y,Z,U,ny,s)
	% Returns the Jacobian matrix (dK/dU)(Y,Z,U), where Y=[Y_1;...;Y_s],
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
	% 	s: number of internal stages
	%
	% Outputs:
	%	dKdU_YZU: The (ny*s)x(s) Jacobian matrix (dK/dU)(Y,Z,U).

	dKdU_YZU = (-1)*blkdiag_vec(Y,ny,s);
end


