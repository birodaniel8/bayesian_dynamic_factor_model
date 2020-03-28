function beta = sampling_beta(x, y, beta_prior, V_prior, sigma_squared, varargin)
%SAMPLING_BETA - Sampling beta from the posterior distribution for normal linear model
%
% Inputs:
%    x: (T x k) explanatory variables
%    y: (T x 1) depentent variabe
%    beta_prior: (double or k x 1) mean of the prior distribution
%    V_prior: (double or k x k) covariance of the prior distribution
%    sigma_sqruared: (double or T x T) specified sigma^2 or covariance matrix of the regression
%    stationary (optional): (bool) sampling betas with stationary constraints
%    last_truncated (optional): (bool) left truncate the last values of the sampled beta at the truncation values
%    truncation_value (optional): (double) applied truncation values of the last element
%
% Outputs:
%    beta: (k x 1) sampled beta values
%
% Note:
%    If beta_prior is a single number a (k x 1) vector is created filled with the given number
%    If V_prior is a single number a (k x k) matrix is created filled with the given number in the diagonals and 0 elsewhere

%------------- BEGIN CODE --------------

%--- INPUT PARSER ---
p = inputParser;
addOptional(p, 'stationary', false, @(x)islogical(x));
addOptional(p, 'last_truncated', false, @(x)islogical(x));
addOptional(p, 'truncation_value', 0, @(x)isnumeric(x));
parse(p, varargin{:});
stationary = p.Results.stationary;
last_truncated = p.Results.last_truncated;
truncation_value = p.Results.truncation_value;
%--------------------

k = size(x, 2);
if all(size(beta_prior) == [1, 1]), beta_prior = repmat(beta_prior, k, 1); end % convert beta_prior to (k x 1)
if all(size(V_prior) == [1, 1]), V_prior = eye(k) * V_prior; end % convert V_prior to (k x k)

V_posterior = inv(inv(V_prior)+(x' * inv(sigma_squared) * x)); % posterior V
beta_posterior = V_posterior * (inv(V_prior) * beta_prior + x' * inv(sigma_squared) * y); % posterior beta

if ~stationary
    if ~last_truncated
        beta = mvnrnd(beta_posterior, V_posterior)';
    else
        if k == 1
            beta = normt_rnd(beta_posterior(k), V_posterior(k, k), truncation_value, Inf);
        else
            beta = [mvnrnd(beta_posterior(1:k-1), V_posterior(1:k-1, 1:k-1)), ...
                normt_rnd(beta_posterior(k), V_posterior(k, k), truncation_value, Inf)]';
        end
    end
else
    stac = 0;
    while stac == 0
        if ~last_truncated
            beta = mvnrnd(beta_posterior, V_posterior)';
        else
            if k == 1
                beta = normt_rnd(beta_posterior(k), V_posterior(k, k), truncation_value, Inf);
            else
                beta = [mvnrnd(beta_posterior(1:k-1), V_posterior(1:k-1, 1:k-1)), ...
                    normt_rnd(beta_posterior(k), V_posterior(k, k), truncation_value, Inf)]';
            end
        end
        new_lagp = LagOp([1, -beta']);
        stac = isStable(new_lagp); % check whether it is stationary
    end
end

if beta(k) == Inf
    disp(beta_posterior)
    disp(V_posterior)
    error("Sampling from the truncated normal distribution was not successful");
end
    %------------- END OF CODE --------------
end