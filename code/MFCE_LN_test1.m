% Test case to find the probability that Z >= t where Z is standard normal.
% P_p(Z >= t) = P_p(f^{(L)}(Z) >= t)

N = 100;
m = 1.6;
s = 0.25;

Z = exp(m + s*randn(N,1));

p = @(x) normpdf(log(x), m, s)./x;

t = 25;
rho = 0.9;
dt = 0.01;

L = 10;
f = @(x, l) x + (exp(-l) - exp(-L));

lvl = 1:L;


[prob, mu, S] = MFCE_LogNormal(Z, t, p, f, lvl, N, rho, dt);

fprintf("\nMFCE Estimate      = %f\n",prob);

prob_true = 1 - normcdf((log(t) - m) / s);
fprintf("True Probability = %f\n", prob_true);

relerr = abs(prob_true - prob)/prob_true;
fprintf("Relative Error   = %f\n\n", relerr);

figure(1);
x = linspace(0.5, 8, 250);
pdf = p(x);
plot(x, pdf, 'b-');
hold on
bias = normpdf(log(x), mu, sqrt(S))./x;
plot(x, bias, 'r-');
legend(["Target", "Biasing"], 'interpreter', 'latex', 'Location', 'NorthWest');
xlabel('$x$', 'interpreter', 'latex');
ylabel('Density', 'interpreter', 'latex');
title('MFCE Biasing Density', 'interpreter', 'latex');