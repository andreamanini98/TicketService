function F = Unif_cdf(x, p)
	a = p(1);
	b = p(2);
	
    % The products here are used as a sort of if..then
	F = max(0, min(1, (x>a) .* (x<b) .* (x - a) / (b - a) + (x >= b)));
end