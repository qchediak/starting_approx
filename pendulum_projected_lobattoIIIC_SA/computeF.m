function F_YZ = computeF(Y,Z,f,n,s)
	% returns F(Y,Z), where Y=[Y_1;...;Y_s], Z=[Z_1;...;Z_s], and
	% F(Y,Z) = [f(Y_1,Z_1);...;f(Y_s,Z_s)].
	%
	% Inputs:
	%	Y: Column vector of internal stages [Y1;...;Ys]
	%	Z: Column vector of internal stages [Z1;...;Zs]
	%	f: RHS of y'=f(y,z)
	%	n: common dimension of Y1,...,Ys, Z1,...,Zs
	% 	s: number of internal stages
	%
	% Outputs:
	%	F_YZ: This is the column vector 
	% 		F(Y,Z) = [f(Y_1,Z_1);...;f(Y_s,Z_s)].

    Y_matrix = reshape(Y, n, s); % convert Y to a matrix [Y1 Y2 ... Ys]
    Z_matrix = reshape(Z, n, s); % convert Z to a matrix [Z1 Z2 ... Zs]
    Y_cells = mat2cell(Y_matrix, n, ones(1, s)); % convert matrix to cell array structure
    Z_cells = mat2cell(Z_matrix, n, ones(1, s)); % convert matrix to cell array structure
    F_cells = cellfun(f, Y_cells, Z_cells, 'UniformOutput', false); % apply f to each pair of cells
    F_YZ = vertcat(F_cells{:}); % concatenate into a single column
end

