function lambda = sampling_mixture_scale_parameter(x, sigma_squared, v)
%SAMPLING_MIXTURE_SCALE_PARAMETER - Sampling the mixture scale parameter
%for the normal linear model with known heteroscedasticity (ie t-errors)
%
% Inputs:
%    x: (T x m) error vector of the regression (conditioned on beta)
%    sigma_squared: (m x 1) static variance of the error term
%    v: (m x 1) degree of freedom parameter
%
% Outputs:
%    lambda: (T x m) sampled mixture scale parameter
%

%------------- BEGIN CODE --------------

[T, m] = size(x);
lambda = zeros(T, m);
for i=1:m
    for t=1:T
        alpha = (v(i) + 1) / 2;
        beta = 2 / ((x(t,i)^2)/sigma_squared(i)+v(i));
        lambda(t,i) = gamrnd(alpha, beta);
    end
end
%------------- END OF CODE --------------
end