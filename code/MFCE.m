function [prob, mu, S] = MFCE(Z, t, p, f, lvl, N, rho, dt, Smin)
% MFCE
%
% Use the MFCE method to estimate the rare-event probability P_p[f(Z) > t].  
% Finds the optimal Gaussian biasing density.
%
% Authors: Terrence Alsup and Frederick Law
%--------------------------------------------------------------------------
% Input:
% 1. Z    -- N-by-d matrix of samples from p.  d is the dimension.
% 2. t    -- the rare-event threshold.
% 3. p    -- the target density as a function handle from R^d to R.
% 4. f    -- the forward model as a function handle from (R^d, lvl) to R.
% 5. lvl  -- the ordered list of levels we evaluate at.
% 6. N    -- the number of samples at each iteration.
% 7. rho  -- the quantile to use ~ 0.9 or 0.99.
% 8. dt   -- the minimum threshold increase at each iteration.
% 9. Smin -- minimum absolute value allowed for the covariance matrix.
%--------------------------------------------------------------------------
% Output:
% 1. prob -- an estimate of the rare-event probability.
% 2. mu   -- the mean of the optimal Guassian biasing density d-by-1.
% 3. S    -- the covariance matrix of the optimal Gaussian biasing density.
%            has size d-by-d.
%--------------------------------------------------------------------------

d = size(Z, 2); % Get the dimension.

% Loop over all levels.
for l = 1:length(lvl)
    mlvl = {}; % Store the parameters from the previous level here.
    if l ~= 1              
        %S = S + Smin*eye(d);
        Z = mvnrnd(mu, S, N);
        mlvl{1} = mu;
        mlvl{2} = S;
    end
    fl = @(x) f(x, lvl(l)); % Model at level l.
    [prob, mu, S] = CrossEntropy(Z, t, p, fl, N, rho, dt, mlvl);
end

end