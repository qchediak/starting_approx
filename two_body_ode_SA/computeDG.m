function J = computeDG(y, G_params)
	% function J = computeDG(y, G_params)
	% Compute the Jacobian DG(y) of the function G(y), where G is the 
	% RHS of the two-body problem y'=G(y).
	% 
	% Inputs: 
	%	y: y=(q,p), where q=(q1,q2) and p=(p1,p2)
	%	G_params: G_params={d}, where d is the common dimension of q and p
	%
	% Outputs: 
	%	J: This is DG(y).

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

	% compute DG as a 2x2 block matrix
	J11 = zeros(d);
	J12 = eye(d);
	J21 = 3*(dot(q,q)^(-5/2))*(q*q') - (dot(q,q)^(-3/2))*eye(d);
	J22 = zeros(d);

	J = [J11 J12; J21 J22];

end

