function F = Erl_cdf(x, p)
	k = p(1);
	lambda = p(2);
	
    % Gammainc computes the regularized lower incomplete gamma function P.
    % The Erlang CDF have been adapted from wikipedia.
    % However, wikipedia seems to have the arguments of P swapped with respect to the
    % ones used by matlab: hence our final formula is as below.
    % To have a quick proof of this claim, simply look at the Gamma argument 
    % and the interval of integration (especially the upper of this interval).
    F = max(0, gammainc(lambda * x, k));
end