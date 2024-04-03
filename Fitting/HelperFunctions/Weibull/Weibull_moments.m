function F = Weibull_moments(p)
	% p(1) -> lambda
	% p(2) -> k
	l = p(1);
	k = p(2);
	
	F = [];
	F(1) = l * gamma(1 + 1 / k);
	F(2) = l^2 * gamma(1 + 2 / k);
end