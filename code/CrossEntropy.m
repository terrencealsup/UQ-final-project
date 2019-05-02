function [prob, mu, S] = CrossEntropy(Xi, rho, N, t, p, F)
% CrossEntropy
%--------------------------------------------------------------------------
% Input:
% 1. Xi  -- a vector of size N of samples from p.
% 2. rho -- the quantile to use.
% 3. N   -- the number of samples to draw at each iteration.
% 4. t   -- the rare-event threshold.
% 5. p   -- the target density.
% 6. F   -- the forward model.
%--------------------------------------------------------------------------
% Output:
% 1. p   -- an estimate of the rare-event probability.
% 2. mu  -- the mean of the optimal Guassian biasing density.
% 3. S   -- the covariance matrix of the optimal Gaussian biasing density.
%--------------------------------------------------------------------------

% The number of samples are rare-events.
Ne = ceil((1 - rho) * N);

Z = Xi;

d = size(Xi, 2); % The dimension of the random variables.

Z = Xi;

k = 0;
while true
    if k == 0
        Fx = F(Xi); % Column vector.
        Gk = sort(Fx);
        gamma = Gk(Ne);
        tk = min([gamma, t]);
        % Fit a Gaussian to the samples.
        indicator = (Fx > tk);
        mu = sum(Xi.*indicator)/sum(indicator);
        S = zeros(d); % d-by-d matrix
        for i=1:N
            if indicator(i)
                S = S + (Xi(i,:)' - mu) * (Xi(i,:)' - mu)';
            end
        end
        S = S / sum(indicator); % Covariance matrix.
        
    else
        Z = mvnrnd(mu, S, N);  % Draw new samples from normal dist.
        Fx = F(Xi);
        Gk = sort(Fx);
        gamma = Gk(Ne);
        tk = min([gamma, t]);
        indicator = (Fx > tk);
        W = p(Z)./mvnpdf(Z, mu, S); % Weights of the samples.
        mu = sum(W'*indicator.*Xi)/sum(W'*indicator);
        S = zeros(d, d);
        for i=1:N
            if indicator(i)
                S = S + (Z(i,:) - mu) * (Z(i,:) - mu)';
            end
        end
        S = S / sum(indicator); % Covariance matrix.
    end

    if tk >= t
        break;
    end
    k = k + 1;
end

prob = mean(W.*indicator);

end

