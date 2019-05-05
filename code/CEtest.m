% Test case to find the probability that Z > t where Z is standard normal.

N = 1e5;
Z = randn(N,1);

p = @(x) 1/sqrt(2*pi) * exp(-0.5*x.^2);
f = @(x) x;

t = 2;
rho = 0.9;
dt = 0.01;

[prob, mu, S] = CrossEntropy(Z, t, p, f, N, rho, dt, {});



prob

prob_true = 1 - normcdf(t)

relerr = abs(prob_true - prob)/prob_true


figure(1);
x = linspace(-4, 4, 500);
pdf = normpdf(x, 0, 1);
plot(x, pdf, 'b-');
hold on
bias = normpdf(x, mu, sqrt(S));
plot(x, bias, 'r-')
