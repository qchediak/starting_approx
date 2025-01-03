function DF_Y = computeDF(Y,Df,n,s)
	% Returns DF(Y), where DF(Y) is the Jacobian of F evaluated at Y.
	%
	% Inputs: 
	%	Y: Column vector of internal stages [Y1;...;Ys]
	%	Df: Jacobian of f, where f is the RHS of y'=f(y)
	%	n: common dimension of Y1,...,Ys
	% 	s: number of internal stages
	%
	% Outputs:
	%	DF_Y: This is DF(Y) = blkdiag(Df(Y1),...,Df(Ys)).

    Y_matrix = reshape(Y, n, s); % convert Y to a matrix [Y1 Y2 ... Ys]
    Y_cells = mat2cell(Y_matrix, n, ones(1, s)); % convert matrix to cell array structure
    DF_cells = cellfun(Df, Y_cells, 'UniformOutput', false); % apply Df to each cell
    DF_Y = blkdiag(DF_cells{:}); % make a block diagonal matrix
end

