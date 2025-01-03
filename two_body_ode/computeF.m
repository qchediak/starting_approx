function F_Y = computeF(Y,f,n,s)
	% Given f, returns F(Y), where Y=[Y_1;...;Y_s] and
	% F(Y) = [f(Y_1);...;f(Y_s)].
	%
	% Inputs:
	%	Y: Column vector of internal stages [Y1;...;Ys]
	%	f: RHS of y'=f(y)
	%	n: common dimension of Y1,...,Ys
	% 	s: number of internal stages
	%
	% Outputs:
	%	F_Y: This is the column vector F(Y) = [f(Y1);...;f(Ys)].

    Y_matrix = reshape(Y, n, s); % convert Y to a matrix [Y1 Y2 ... Ys]
    Y_cells = mat2cell(Y_matrix, n, ones(1, s)); % convert matrix to cell array structure
    F_cells = cellfun(f, Y_cells, 'UniformOutput', false); % apply f to each cell
    F_Y = vertcat(F_cells{:}); % concatenate into a single column
end

