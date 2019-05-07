addpath('../ellip/');


N = 1e4;
t = 0.26;

d = 16;



lvl = 3:5;
rho = 0.9;
dt = 0.01;

m = 1.6 * ones(1,d);
s = 0.25^2 * eye(d);

Z = exp(mvnrnd(m, s, N));


p = @(x) temp_p(x, m, s);

[prob, mu, S] = MFCE_LogNormal(Z, t, p, @f, lvl, N, rho, dt);

prob
true_prob = 5.4e-3;

relerr = abs(true_prob - prob)/true_prob;

function vals = f(x, l)
    ops = load(['operatorsBlocks16_level',num2str(l),'.mat']);
    vals = ellip2DAffine(x,ops.ACell,ops.f,l);
end

function vals = temp_p(x,m,s)
    vals = mvnpdf(log(x), m, s) ./ prod(x, 2);
end