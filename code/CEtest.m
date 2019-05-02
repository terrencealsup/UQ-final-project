N = 100;
Xi = mvnrnd(0,1,N);

p = @(x) 1/sqrt(2*pi) * exp(-0.5*x.^2);
F = @(x) x;

CrossEntropy(Xi, 0.1, N, 2, p, F)