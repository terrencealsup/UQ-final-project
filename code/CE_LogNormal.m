function [prob, mu, S] = CE_LogNormal(Z, t, p, f, N, rho, dt, mlvl)
% CrossEntropy
%
% Use the Cross Entropy method to estimate the rare-event probability
% P_p[f(Z) > t].  Finds the optimal log-normal biasing density.
%
% Authors: Terrence Alsup and Frederick Law
%--------------------------------------------------------------------------
% Input:
% 1. Z    -- N-by-d matrix of samples from p.  d is the dimension.
% 2. t    -- the rare-event threshold.
% 3. p    -- the target density as a function handle from R^d to R.
% 4. f    -- the forward model as a function handle from R^d to R.
% 5. N    -- the number of samples at each iteration.
% 6. rho  -- the quantile to use ~ 0.9 or 0.99.
% 7. dt   -- the minimum threshold increase at each iteration.
% 8. mlvl -- Is this function called from a multilevel setting?
%            If not then mlvl is an empty cell, 
%            otherwise contains the mean mu in first cell as a column
%            vector and the covariance matrix S in the second cell from the
%            previous level.
%--------------------------------------------------------------------------
% Output:
% 1. prob -- an estimate of the rare-event probability.
% 2. mu   -- mean of the optimal Log-Normal biasing density d-by-1.
% 3. S    -- covariance matrix of the optimal Log-Normal biasing density.
%            has size d-by-d.
%--------------------------------------------------------------------------


d  = size(Z, 2); % The dimension of the random variables.

% Initialize relevant variables.
k  = 1;             % Iteration counter.
tk = t - 1;         % Initialize the threshold to be less than t at first.

if isempty(mlvl)
    mu = zeros(d, 1);   % Initialize the mean.
    S  = eye(d);        % Initialize the covariance matrix.
else
    mu = mlvl{1};       % Optimal mean from previous level.
    S  = mlvl{2};       % Optimal covariance from previous level.
end

P  = zeros(N, 1);   % The target density evaluations.
Q  = ones(N, 1);    % The biasing density evaluations.
buffer = 1e-3 * eye(d);

% The main iteration.
while tk < t
    
    % Draw samples based on previous iteration.
    if k ~= 1
        Z = exp(mvnrnd(mu, S, N));  % Draw samples from log-normal.
    end
    
    % Compute model outputs.  Stored as a column vector N-by-1.
    F = f(Z);
    
    % Estimate the rho quantile.
    gamma = quantile(F, rho);
    
    % Update the threshold.
    if k == 1
        tk = min([gamma, t]);
    else  
        tk = min([t, max([tk + dt, gamma])]);
    end

    I = (F >= tk); % The indicator function for each sample.
    
    % Now compute the weights.
    W = ones(N, 1);
    if k ~= 1 || ~isempty(mlvl)
        % Compute the target and biasing densities at the samples.
%         for i = 1:N
%             if I(i) ~= 0
%                 P(i) = p(Z(i,:));
%                 Q(i) = mvnpdf(log(Z(i,:)'), mu, S) ./ prod(Z(i,:));
%                 W(i) = P(i) / Q(i);
%             else
%                 W(i) = 0;
%             end
%         end
        %vectorize
        P = p(Z);
        Q = mvnpdf(log(Z),mu',S) ./ prod(Z,2);
        W = P ./ Q;
    end
    
    % Estimate the mean.
%     mu = zeros(d, 1);
%     for i = 1:N
%         if I(i) ~= 0 && W(i) ~= 0
%             mu = mu + (I(i) * W(i) * log(Z(i,:)'));
%         end
%     end
%     mu = mu / (I' * W); % Re-normalize.

    % Vectorized.
    temp1 = I .* W; 
    sum1 = sum(temp1);
    mu = (log(Z)' * temp1) / sum1;
    
    % Estimate the covariance matrix.
    S = zeros(d);
    S_temp = log(Z)' - mu*ones(1,N);
    for i = 1:N
        if I(i) ~= 0 && W(i) ~= 0
%             S = S + I(i) * W(i) * ((log(Z(i,:)') - mu) * (log(Z(i,:)') - mu)');
            temp2 = S_temp(:,i);
            S = S + temp1(i) * (temp2 * temp2');
        end
    end
%     S = S / (I' * W); % Re-normalize.
    S = S / sum1;
    
    % Inflate the covariance matrix so it is better conditioned.
%     [V, D] = eig(S);
%     D = max(D, 1e-3 * eye(d)); % Minimum eigenvalue is at least 10^-3.
%     S = V * D * V';
%     S = S + 1e-3 * eye(d);
      S = S + buffer;
    
    % Update iteration counter.
    k = k + 1;
end
disp(['# CE iterations is K=',num2str(k-1)]);
% Now compute the Cross Entropy estimate for the probability.

% prob = (I' * W) / N;
prob = sum1 / N;

end
