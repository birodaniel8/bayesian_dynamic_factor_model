function v = sampling_sigma_squared(x, gamma_prior, delta_prior)
%SAMPLING_SIGMA_SQUARED - Sampling sigma^2 from the posterior distribution for normal linear model
%
% Inputs:
%    x: (T x m) error vector of the regression (conditioned on beta)
%    gamma_prior: (m x 1) shape of the prior distribution
%    delta_prior: (m x 1) scale of the prior distribution
%
% Outputs:
%    v: (m x 1) sampled sigma^2 values
%

%------------- BEGIN CODE --------------

% we consider only independent sigma squared priors, therefor we take only the diagonal of the x'*x matrix
v_posterior = size(x,1)/2 + gamma_prior;
d_posterior = diag(x'*x)/2 + delta_prior;
v = 1 ./ gamrnd(v_posterior, 1./d_posterior);

%------------- END OF CODE --------------   
end
