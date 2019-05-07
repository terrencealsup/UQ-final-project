addpath('../ellip/');

p = @(x) ones(size(x)).*(x >= 1).*(x<=10)/9;

N = 1e3;
t = 0.9;

Z = 9*rand(N, 1) + 1;

lvl = 3:5;
rho = 0.9;
dt = 0.01;


[prob, mu, S] = MFCE(Z, t, p, @f, lvl, N, rho, dt);

prob


function vals = f(x, l)
    ops = load(['operatorsBlocks1_level',num2str(l),'.mat']);
    vals = ellip2DAffine(x,ops.ACell,ops.f,l);
end