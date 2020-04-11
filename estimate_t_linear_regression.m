function [sampled_beta, sampled_sigma_squared, sampled_df] = estimate_t_linear_regression(y,x,priors,ndraw,burnin,constant,display)
% ESTIMATE_T_LINEAR_REGRESSION

    if isfield(priors,'beta'), beta_prior = priors.beta; else beta_prior = 0; end  % prior mean of the beta coefficients in the normal model
    if isfield(priors,'V'), V_prior = priors.V; else V_prior = 1; end  % prior variance of the beta coefficients in the normal model
    if isfield(priors,'gamma'), gamma_prior = priors.gamma; else gamma_prior = 3/2; end  % prior confidence in the variance in the normal model
    if isfield(priors,'delta'), delta_prior = priors.delta; else delta_prior = 0.01/2; end  % prior mean of the variance in the normal model
    if isfield(priors,'v'), v_prior = priors.v; else v_prior = 30; end  % prior mean of the degree of freedom parameter
    if isfield(priors,'mh_variance'), mh_variance = priors.mh_variance; else mh_variance = 0.25; end  % prior mean of the variance in the normal model

    [T, m] = size(x);
    if constant
        x = [ones(T,1) x];
        m = m + 1;
    end
    
    % create containers:
    sampled_beta = zeros(m, ndraw-burnin);
    sampled_sigma_squared = zeros(1, ndraw-burnin);
    sampled_df = zeros(1, ndraw-burnin);
    
    % initial values:
    sigma_squared = delta_prior;
    lambda = ones(T,1);
    df = v_prior;

    % sampling:
    for i = 1:ndraw
        if display && mod(i,100) == 0, disp(i); end
        x_star = sqrt(lambda) .* x;
        y_star = sqrt(lambda) .* y;
        beta = sampling_beta(x_star, y_star, beta_prior, V_prior, sigma_squared);   % sampling beta coefficients
        sigma_squared = sampling_sigma_squared(y_star - x_star * beta, gamma_prior, delta_prior);  % sampling error variance
        df = sampling_degree_of_freedom(lambda, df, v_prior, mh_variance);  % sampling degree of freedom
        lambda = sampling_mixture_scale_parameter(y - x * beta, sigma_squared, df);  % sampling lambda
        if i > burnin
            sampled_beta(:, i-burnin) = beta;
            sampled_sigma_squared(:, i-burnin) = sigma_squared;
            sampled_df(:, i-burnin) = df;            
        end
    end
end