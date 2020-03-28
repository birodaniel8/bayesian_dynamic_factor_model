function v = sampling_degree_of_freedom(lambda, v_previous, v_prior, mh_sigma_squared)
%SAMPLING_DEGREE_OF_FREEDOM - Sampling the degree of freedom parameter
%for the normal linear model with known heteroscedasticity (ie t-errors)
%
% Inputs:
%    lambda: (T x m) mixture scale parameters
%    v_previous ((m) x 1) vector of the previous v parameters
%    v_prior: ((m) x 1) vector of prior degree of freedom parameters
%    mh_sigma_squared: ((m) x 1) random walk Metropolis-Hastings algorithm
%                                variance parameters
%
% Outputs:
%    v: (m x 1) sampled degree of freedom parameter
%

%------------- BEGIN CODE --------------

m = size(lambda,2);
if all(size(v_previous) == [1, 1]), v_previous = repmat(v_previous, m, 1); end % convert v_previous to (m x 1)
if all(size(v_prior) == [1, 1]), v_prior = repmat(v_prior, m, 1); end % convert v_prior to (m x 1)
if all(size(mh_sigma_squared) == [1, 1]), mh_sigma_squared = repmat(mh_sigma_squared, m, 1); end % convert mh_sigma_squared to (m x 1)

v = zeros(m, 1);
for i=1:m
    v_proposed = v_previous(i) + sqrt(mh_sigma_squared(i)) * randn(1);
    alpha = min(v_posterior_density(v_proposed,v_prior(i),lambda(:,i))/v_posterior_density(v_previous(i),v_prior(i),lambda(:,i)),1);
    if rand() < alpha
        v(i) = v_proposed;
    else
        v(i) = v_previous(i);
    end
end

%------------- END OF CODE --------------
end

%------- BEGIN AUXILIARY FUNCTION -------
function x = v_posterior_density(v,v_prior,lambda)
%v_POSTERIOR_DENSITY: p(v|v_prior,N,lambda)
    [N,~] = size(lambda);
    eta = 1/v_prior + 1/2 * sum(log(1./lambda) + lambda);
    x = (v/2)^(N*v/2) * gamma(v/2)^(-N) * exp(-eta*v);
end
%-------- END AUXILIARY FUNCTION --------