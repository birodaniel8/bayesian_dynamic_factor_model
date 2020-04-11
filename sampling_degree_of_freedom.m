function v = sampling_degree_of_freedom(lambda, v_previous, v_prior, mh_variance)
%SAMPLING_DEGREE_OF_FREEDOM - Sampling the degree of freedom parameter
%for the normal linear model with known heteroscedasticity (ie t-errors)
%
% Inputs:
%    lambda: (T x m) mixture scale parameters
%    v_previous ((m) x 1) vector of the previous v parameters
%    v_prior: ((m) x 1) vector of prior degree of freedom parameters
%    mh_variance: ((m) x 1) random walk Metropolis-Hastings algorithm
%                           variance parameters
%
% Outputs:
%    v: (m x 1) sampled degree of freedom parameter
%

%------------- BEGIN CODE --------------

[N, m] = size(lambda);
if all(size(v_previous) == [1, 1]), v_previous = repmat(v_previous, m, 1); end % convert v_previous to (m x 1)
if all(size(v_prior) == [1, 1]), v_prior = repmat(v_prior, m, 1); end % convert v_prior to (m x 1)
if all(size(mh_variance) == [1, 1]), mh_variance = repmat(mh_variance, m, 1); end % convert mh_sigma_squared to (m x 1)

v = zeros(m, 1);
for i=1:m
    v_proposed = v_previous(i) + sqrt(mh_variance(i)) * randn(1);

    eta =1/v_prior + 1/2 *sum(-log(lambda(:,i)) + lambda(:,i));
    
    if v_proposed>0
        lpostcan = .5*N*v_proposed*log(.5*v_proposed) -N*gammaln(.5*v_proposed)-eta*v_proposed;
        lpostdraw = .5*N*v_previous(i)*log(.5*v_previous(i)) -N*gammaln(.5*v_previous(i))-eta*v_previous(i);
        accprob = exp(lpostcan-lpostdraw);
     else
        accprob=0;
     end
    
    alpha = accprob;
    if rand < alpha
        v(i) = v_proposed;
    else
        v(i) = v_previous(i);
    end
end

%------------- END OF CODE --------------
end