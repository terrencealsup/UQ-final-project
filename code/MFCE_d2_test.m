addpath('../ellip/');


N = 1e5;
t = 0.45;

%Z = rand(N, 16);
d = 16;
Z = mvnrnd(zeros(d,1), 1e-2*eye(d), N);


lvl = 3:5;
rho = 0.9;
dt = 0.01;


[prob, mu, S] = MFCE(Z, t, @p, @f, lvl, N, rho, dt);


true_prob = 1e-4;

relerr = abs(true_prob - prob)/true_prob;

function vals = f(x, l)
    ops = load(['operatorsBlocks16_level',num2str(l),'.mat']);
    x = 9*x + 1;
    vals = ellip2DAffine(x,ops.ACell,ops.f,l);
end

function vals = p(x)
    vals = (prod(x >= zeros(size(x)) & x <= ones(size(x))) == 1);
end