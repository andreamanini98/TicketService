function F = Pareto_momentsEquation(p)
	global M1_par
	global M2_par
	% p(1) -> alpha
	% p(2) -> m
	
	F = Pareto_moments(p);
	F(1) = F(1) / M1_par - 1; % alpha
	F(2) = F(2) / M2_par - 1; % m
end