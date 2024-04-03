% Change the name in csvread() to select a different trace.
Trace = csvread('TraceC-I.txt');
N = size(Trace,1);

% Since times in traces are expressed in milliseconds, in order to get parameters 
% in seconds, we divide them by 1000 [ 1 millisecond = 0.001 second ].
% If you only want to show plots, comment out the division.
Trace = Trace / 1000;

% The range in which to plot our fitted distribution.
% Change this accordingly, in order to have a decent plot based on the trace you are analyzing.
Range = 200;

% Finding the first three moments.
M1 = sum(Trace) / N;
M2 = sum(Trace .^ 2) / N;
M3 = sum(Trace .^ 3) / N;

% Computing standard deviation and coefficient of variation of the samples.
Sigma = std(Trace);
Cv = Sigma / M1;

SortedTrace = sort(Trace);


% Now fitting Exponential with direct expression.
lambda_exp = 1 / M1;
% t is the range for plotting our Exponential.
t_exp = [1:Range] / 10;
% Now plotting the sampled CDF from the estimated Exponential one.
figure('NumberTitle', 'off', 'Name', 'Exponential direct fitting');
p = plot(SortedTrace, [1:N]/N, "-", t_exp, Exp_cdf(t_exp, [lambda_exp]), "-");
grid on;
legend({'Samples','Fitting'},'Location','southeast');
p(1).LineWidth = 1.5;
p(2).LineWidth = 1.5;


% Now fitting Uniform with direct expression.
a_uni = M1 - sqrt(12 * (M2 - M1^2)) / 2;
b_uni = M1 + sqrt(12 * (M2 - M1^2)) / 2;
% t is the range for plotting our Uniform.
t_uni = [1:Range] / 10;
% Now plotting the sampled CDF from the estimated Uniform one.
figure('NumberTitle', 'off', 'Name', 'Uniform direct fitting');
p = plot(SortedTrace, [1:N]/N, "-", t_uni, Unif_cdf(t_uni, [a_uni, b_uni]), "-");
grid on;
legend({'Samples','Fitting'},'Location','southeast');
p(1).LineWidth = 1.5;
p(2).LineWidth = 1.5;


% Now fitting Erlang with direct expression.
% Since we know that: mu = k / l and var = k / l^2.
% We obtain that: k = mu * l.
% Substituting in the second equation: var = (mu * l) / l^2 = mu / l.
% From this we can obtain: l = mu / var.
% Then, substituting this in the first equation we get: k = mu * (mu / var) = mu^2 / var.
% k must be an integer since by definition an Erlang distribution is a summation of k exponentials.
k_erl = round(M1^2 / Sigma^2);
lambda_erl = M1 / Sigma^2;
% t is the range for plotting our Erlang.
t_erl = [1:Range] / 10;
% Now plotting the sampled CDF from the estimated Erlang one.
figure('NumberTitle', 'off', 'Name', 'Erlang direct fitting');
p = plot(SortedTrace, [1:N]/N, "-", t_erl, Erl_cdf(t_erl, [k_erl, lambda_erl]), "-");
grid on;
legend({'Samples','Fitting'},'Location','southeast');
p(1).LineWidth = 1.5;
p(2).LineWidth = 1.5;


% Now fitting two stage Hypo-Exponential with Maximum Likelyhood.
% What to pass to mle:
% - The dataset
% - 'pdf' keyword
% - A function description of the pdf: @(dataset, parameters_list)Function_name(dataset, [parameters_list])
% - 'start' keyword
% - A vector containing the starting point
par_hypo = mle(Trace, 'pdf', @(Trace, l1, l2)HypoExp_pdf(Trace, [l1, l2]), 'start', [1 / (0.3 * M1), 1 / (0.7 * M1)]);
% t is the range for plotting our Hypo-Exponential.
t_hypo = [1:Range] / 10;
% Now plotting the sampled CDF from the estimated Hypo-Exponential one.
figure('NumberTitle', 'off', 'Name', 'Hypo-Exponential Maximum Likelyhood fitting');
p = plot(SortedTrace, [1:N]/N, "-", t_hypo, HypoExp_cdf(t_hypo, par_hypo), "-");
grid on;
legend({'Samples','Fitting'},'Location','southeast');
p(1).LineWidth = 1.5;
p(2).LineWidth = 1.5;


% Now fitting two stage Hyper-Exponential with Maximum Likelyhood.
% WARNING: THE COEFFICIENT OF VARIATION IS < 1, SO THE ESTIMATION WILL BE WRONG!!
par_hyper = mle(Trace, 'pdf', @(Trace, l1, l2, p1)HyperExp_pdf(Trace, [l1, l2, p1]), 'start', [0.8 / M1, 1.2 / M1, 0.4]);
% t is the range for plotting our Hyper-Exponential.
t_hyper = [1:Range] / 10;
% Now plotting the sampled CDF from the estimated Hyper-Exponential one.
figure('NumberTitle', 'off', 'Name', 'Hyper-Exponential Maximum Likelyhood fitting');
p = plot(SortedTrace, [1:N]/N, "-", t_hyper, HyperExp_cdf(t_hyper, par_hyper), "-");
grid on;
legend({'Samples','Fitting'},'Location','southeast');
p(1).LineWidth = 1.5;
p(2).LineWidth = 1.5;


% Now fitting Weibull with Method of Moments.
% We have to decleare two global variables to be used in Weibull_momentsEquation().
global M1_weib;
global M2_weib;
M1_weib = M1;
M2_weib = M2;
% [1, 1] represent the fitting starting point.
% par_weib(1) = lambda, par_weib(2) = k.
par_weib = fsolve(@(x)Weibull_momentsEquation(x), [1, 1]);
% t is the range for plotting our Weibull.
t_weib = [1:Range] / 10;
% Now plotting the sampled CDF from the estimated Weibull one.
figure('NumberTitle', 'off', 'Name', 'Weibull Method of Moments fitting');
p = plot(SortedTrace, [1:N]/N, "-", t_weib, Weibull_cdf(t_weib, par_weib), "-");
grid on;
legend({'Samples','Fitting'},'Location','southeast');
p(1).LineWidth = 1.5;
p(2).LineWidth = 1.5;


% Now fitting Pareto with Method of Moments.
% We have to decleare two global variables to be used in Pareto_momentsEquation().
global M1_par;
global M2_par;
M1_par = M1;
M2_par = M2;
% [3, 3] represent the fitting starting point. We chose [3, 3] since alpha must be at least > 2
% in order not to have infinty-valued moments.
% par_par(1) = alpha, par_par(2) = m.
par_par = fsolve(@(x)Pareto_momentsEquation(x), [3, 3]);
% t is the range for plotting our Pareto.
t_par = [1:Range] / 10;
% Now plotting the sampled CDF from the estimated Pareto one.
figure('NumberTitle', 'off', 'Name', 'Pareto Method of Moments fitting');
p = plot(SortedTrace, [1:N]/N, "-", t_par, Pareto_cdf(t_par, par_par), "-");
grid on;
legend({'Samples','Fitting'},'Location','southeast');
p(1).LineWidth = 1.5;
p(2).LineWidth = 1.5;
