function F = Pareto_cdf(x, p)
	alpha = p(1);
	m = p(2);

	F = max(0, (x >= m) .* (1 - power(m ./ x, alpha)));
end