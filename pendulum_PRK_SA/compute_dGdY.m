function dGdY_Y = compute_dGdY(Y,ny,s)
	% function dGdY_Y = compute_dGdY(Y,ny,s)
	% Returns the Jacobian matrix (dG/dY)(Y), where Y=[Y_1;...;Y_s], and
	% G(Y)=[g(Y1);...;g(Ys)].
	%
	% This is specific to the pendulum problem.
	%
	% Inputs:
	%	Y: Column vector of internal stages [Y1;...;Ys]
	%	ny: common dimension of Y1,...,Ys
	% 	s: number of internal stages
	%
	% Outputs:
	%	dGdY_Y: The sx(ny*s) matrix (dG/dY)(Y)

	% split Y into a cell array of ny x 1 vectors
	Y_cells = mat2cell(Y, ny * ones(s,1), 1);

	% apply transpose and scaling to each block
	Yt_cells = cellfun(@(y) 2*y', Y_cells, 'UniformOutput', false);

	% Make the block diagonal matrix
	dGdY_Y = blkdiag(Yt_cells{:});
end


