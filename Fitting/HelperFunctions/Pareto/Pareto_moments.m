function F = Pareto_moments(p)
	% p(1) -> alpha
	% p(2) -> m
	a = p(1);
	m = p(2);
	
	F = [];
    if (a <= 1)
        F(1) = Inf;
    else
        F(1) = a * m / (a - 1);
    end

    if (a <= 2)
        F(2) = Inf;
    else
        F(2) = a * m^2 / (a - 2);
    end
end