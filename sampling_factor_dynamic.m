function factor = sampling_factor_dynamic(y, beta, theta, error_covariance, factor_covariance)
%SAMPLING_FACTOR_NORMAL_DYNAMIC - Sampling factor values from the posterior distribution for normal dynamic factor model,
% where the factors follow 
%
% Inputs:
%    y: (T x m) matrix of observations
%    beta: (m x k x (T)) factor loadings
%    theta: (k x 1) autoregressive coefficient of the factors
%    error_covariance: (m x m x (T)) covariance matrix of the observation equation
%    factor_covariance: (k x k x (T)) covariance matrix of the factor equation
%
% Outputs:
%    factor: (T x k) sampled factor values

%------------- BEGIN CODE --------------

k = size(beta, 2);
if nargin == 4, factor_covariance = eye(k); end

H = beta;
R = error_covariance;
G = diag(theta);
Q = factor_covariance;
x0 = zeros(k,1);
P0 = eye(k);

[F, P] = kalman_filter(y, H, R, G, Q, x0, P0);
factor = sampling_carter_kohn(F, P, G, Q);

%------------- END OF CODE --------------
end