function factor = sampling_factor_normal(y, beta, error_covariance)
%SAMPLING_FACTOR_NORMAL - Sampling factor values from the posterior distribution for normal linear model (f ~ N(0,1))
%
% Inputs:
%    y: (T x m) matrix of observations
%    beta: (m x k) factor loadings
%    error_covariance: (m x m) covariance matrix of the observation equation
%
% Outputs:
%    factor: (T x K) sampled factor values

%------------- BEGIN CODE --------------

T = size(y, 1);
k = size(beta, 2);
factor = zeros(T, k);

for t = 1:T
    factor_var = inv(eye(k)+beta'*inv(error_covariance)*beta);
    factor_mean = factor_var * beta' * inv(error_covariance) * y(t, :)';
    factor_sampled = factor_mean + chol(factor_var)' * randn(k, 1);
    factor(t, 1:k) = factor_sampled';
end

%------------- END OF CODE --------------
end