% Test case to find the probability that Z >= t where Z is standard normal.
% P_p(Z >= t) = P_p(f^{(10)}(Z) >= t)

N = 1e5;
Z = randn(N,1);

L = 3; % High-fidelity.

p = @(x) 1/sqrt(2*pi) * exp(-0.5*x.^2);
f = @(x, l) x + (exp(-l) - exp(-L));


lvl = 1:L;
t = 2;
rho = 0.9;
dt = 0.01;
Smin = 1e-3;

[prob, mu, S] = MFCE(Z, t, p, f, lvl, N, rho, dt, Smin);

prob

prob_true = 1 - normcdf(t)

relerr = abs(prob_true - prob)/prob_true