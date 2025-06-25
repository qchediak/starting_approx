function B = blkdiag_vec(V, n, s)
	% function B = blkdiag_vec(V, n, s)
	%
	% Utility function for converting the vector V=[V1;...;Vs]
	% to a block diagonal matrix blkdiag(V1,...,Vs).
	%
	% Inputs:
	%	V: Column vector of [V1;...;Vs] \in \mathbb{R}^{ns}, where each
	%		component V_i is itself a column vector in \mathbb{R}^n
	%	n: The common dimension of V1,...,Vs
	% 	s: The number of components V1,...,Vs
	%
	% Output:
	%	B: The (n*s)xs block diagonal matrix blkdiag(V1,...,Vs)

	U_cells = mat2cell(V, n * ones(s,1), 1);
	B = blkdiag(U_cells{:});
end

