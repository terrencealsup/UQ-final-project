% Plot the biasing density as well as the optimal biasing density for the
% 1-dimensional problem.

addpath('../ellip/');

d = 1; % The dimension.
L = 5; % The level to solve the PDE at.

% Parameters for the Cross Entropy method.
N = 1e4;    % Number of samples.
t = 0.36;   % Rare-event threshold.
rho = 0.9;
dt = 0.01;

% Parameters of the log-normal distribution.
m = 1.6;        % Mean.
s = 0.0625;     % Variance.

% Generate samples from the target density.
Z = exp(m + sqrt(s)*randn(N,1));

% Evaluate the log-normal density.
p = @(x) temp_p(x, m, s);
f = @(x) temp_f(x, L);

% Run the Cross Entropy method.
[prob, mu, S] = CE_LogNormal(Z, t, p, f, N, rho, dt, {});


x = linspace(1, 9, 250)';

% The computed biasing density.
bias = temp_p(x, mu, S);

% The optimal biasing density.
opt_bias = p(x) .* (f(x) >= t) / prob;

% Plot the densities.
figure(1);
plot(x, p(x), 'm-', 'Linewidth', 2);
hold on
plot(x, bias, 'b-', 'Linewidth', 2);
plot(x, opt_bias, 'r-', 'Linewidth', 2);
xlabel('$\xi$', 'interpreter', 'latex');
ylabel('Density', 'interpreter', 'latex');
title('Cross-Entropy and Optimal Biasing Densities', 'interpreter', 'latex');
legend(["Target", "CE", "Optimal"], 'interpreter', 'latex');
hold off



function vals = temp_f(x, l)
    ops = load(['operatorsBlocks1_level',num2str(l),'.mat']);
    vals = ellip2DAffine(x,ops.ACell,ops.f,l);
end

function vals = temp_p(x,m,s)
    vals = normpdf(log(x), m, sqrt(s)) ./ x;
end