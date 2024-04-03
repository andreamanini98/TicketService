function F = Weibull_momentsEquation(p)
	global M1_weib
	global M2_weib
	% p(1) -> lambda
	% p(2) -> k
	
	F = Weibull_moments(p);
	F(1) = F(1) / M1_weib - 1; % lambda
	F(2) = F(2) / M2_weib - 1; % k
end