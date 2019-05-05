% Test case to find the probability that Z >= t where Z is standard normal.
% P_p(Z >= t) = P_p(f^{(L)}(Z) >= t)

N = 1e3;
Z = randn(N,1);

L = 10; % High-fidelity.

p = @(x) 1/sqrt(2*pi) * exp(-0.5*x.^2);
f = @(x, l) x + (exp(-l) - exp(-L));


lvl = 1:L;
t = 2;
rho = 0.9;
dt = 0.01;

[prob, mu, S] = MFCE(Z, t, p, f, lvl, N, rho, dt);

fprintf("\nMFCE Estimate    = %f\n",prob);

prob_true = 1 - normcdf(t);
fprintf("True Probability = %f\n", prob_true);

relerr = abs(prob_true - prob)/prob_true;
fprintf("Relative Error   = %f\n\n", relerr);

figure(1);
x = linspace(-4, 4, 500);
pdf = normpdf(x, 0, 1);
plot(x, pdf, 'b-');
hold on
bias = normpdf(x, mu, sqrt(S));
plot(x, bias, 'r-');
legend(["Target", "Biasing"], 'interpreter', 'latex', 'Location', 'NorthWest');
xlabel('$x$', 'interpreter', 'latex');
ylabel('Density', 'interpreter', 'latex');
title('MFCE Biasing Density', 'interpreter', 'latex');
