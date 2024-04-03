function F = Weibull_cdf(x, p)
	l = p(1);
	k = p(2);

	F = max(0, 1 - exp(-power(x / l, k)));
end