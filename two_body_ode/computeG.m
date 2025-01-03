function G_y = computeG(y,G_params)
	% function G_y = G1(y,G_params)
	% Find G(y), where G is the RHS of the two-body system y'=G(y), 
	% 
	% Inputs: 
	%	y: y=(q,p), where q=(q1,q2) and p=(p1,p2)
	%	G_params: G_params={d}, where d is the common dimension of q and p
	%
	% Outputs: 
	%	G_y: This is G(y) given the value d.

	% unpack G_params
	d = G_params{1};

	% check dimensions of y
	if size(y,1) ~= 2*d
		error('Error (G1): y has the wrong number of rows')
	elseif size(y,2) > 1
		error('Error (G1): y is not a column vector')
	end

	% unpack y
	q = y(1:d);
	p = y(d+1:end);

	% compute G
	G_y = [p; -(dot(q,q)^(-3/2))*q];

end


