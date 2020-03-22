function beta = sampling_factor_loading_normal(y, factor, beta_prior, V_prior, error_covariance)
%SAMPLING_FACTOR_LOADING_NORMAL - Sampling factor loadings from the posterior distribution for normal linear model
%
% Inputs:
%    y: (T x m) matrix of observations
%    factor: (T x k) matrix of factors
%    beta_prior: (double or m x k) mean of the prior distribution
%    V_prior: (double or m x k) variance of the prior distribution
%    error_covariance: (m x m) covariance matrix of the observation equation
%
% Outputs:
%    beta: (m x k) sampled factor loadings (lower triangle)
%
% Note:
%    Factor loadings are sampled from independent normal distributions using N(beta_prior,V_prior) prior
%    If beta_prior is a single number a (m x k) matrix is created filled with the given number
%    If V_prior is a single number a (m x k) matrix is created filled with the given number

%------------- BEGIN CODE --------------

m = size(y, 2);
k = size(factor, 2);
beta = zeros(m, k);

if all(size(beta_prior) == [1, 1]), beta_prior = repmat(beta_prior, m, k); end % convert beta_prior to (k x 1)
if all(size(V_prior) == [1, 1]), V_prior = repmat(V_prior, m, k); end % convert V_prior to (k x k)


for i = 1:m
    if i <= k
        beta(i, 1:i) = sampling_beta(factor(:, 1:i), y(:, i), beta_prior(i, 1:i)', diag(V_prior(i, 1:i)), error_covariance(i, i), 'last_truncated', true);
    else
        beta(i, :) = sampling_beta(factor, y(:, i), beta_prior(i, :)', diag(V_prior(i, :)), error_covariance(i, i));
    end
end

%------------- END OF CODE --------------
end