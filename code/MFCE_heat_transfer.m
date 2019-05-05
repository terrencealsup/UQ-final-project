addpath('../ellip/');

global d = 2;


N = 1e3;
t = 0.78;

Z = 9*rand(N, d) + 1;

lvl = 3:5;
rho = 0.9;
dt = 0.01;


[prob, mu, S] = MFCE(Z, t, @p, @f, lvl, N, rho, dt);



function vals = f(x, l)
    ops = load(['operatorsBlocks',num2str(d),'_level',num2str(l),'.mat']);
    vals = ellip2DAffine(x,ops.ACell,ops.f,l);
end

function vals = p(x)
    vals = ((1/9)^d)*(prod(x >= ones(size(x)) & x <= 10*ones(size(x))) == 1);
end