% Test case to find the probability that Z > t where Z is standard normal.

N = 1e7;
Z = randn(N,1);

p = @(x) 1/sqrt(2*pi) * exp(-0.5*x.^2);
f = @(x) x;

t = 4;
rho = 0.9;
dt = 0.01;

[prob, mu, S] = CrossEntropy(Z, t, p, f, N, rho, dt);

prob_true = 1 - normcdf(t)

relerr = abs(prob_true - prob)/prob_true