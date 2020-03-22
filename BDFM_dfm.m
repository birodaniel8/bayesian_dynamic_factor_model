classdef BDFM_dfm < BDFM
    %BDFM_dfm Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        beta_prior, V_prior, gamma_prior, delta_prior, theta_prior, theta_V_prior
        sampled_beta, sampled_sigma_squared, sampled_factor, sampled_theta
    end
    
    methods
        function obj = BDFM_dfm(y, k, priors)
            %BDFM_nlm Construct an instance of this class
            obj@BDFM(y, k);
            if isfield(priors,'beta'), obj.beta_prior = priors.beta; else obj.beta_prior = 0; end  % prior mean of the beta coefficients in the normal model
            if isfield(priors,'V'), obj.V_prior = priors.V; else obj.V_prior = 1; end  % prior variance of the beta coefficients in the normal model
            if isfield(priors,'gamma'), obj.gamma_prior = priors.gamma; else obj.gamma_prior = 3/2; end  % prior confidence in the variance in the normal model
            if isfield(priors,'delta'), obj.delta_prior = priors.delta; else obj.delta_prior = 0.01/2; end  % prior mean of the variance in the normal model
            if isfield(priors,'theta'), obj.theta_prior = priors.theta; else obj.theta_prior = 0; end  % prior mean of the autogregressive coefficient of the factor equation
            if isfield(priors,'theta_V'), obj.theta_V_prior = priors.theta_V; else obj.theta_V_prior = 1; end  % prior variance of the autogregressive coefficient of the factor equation
        end
        
        function obj = estimate(obj, ndraw, burnin, display)
            if nargin == 3, display = true; end
            % create containers:
            obj.sampled_beta = zeros(obj.m, obj.k, ndraw-burnin);
            obj.sampled_sigma_squared = zeros(obj.m, ndraw-burnin);
            obj.sampled_theta = zeros(obj.k, ndraw-burnin);
            obj.sampled_factor = zeros(obj.T, obj.k, ndraw-burnin);
            % initial values:
            sigma_squared = repmat(obj.delta_prior, obj.m, 1);
            factor = factor_initialize(obj.y,obj.k);
            theta = repmat(obj.theta_prior, obj.k, 1);
            % sampling:
            for i = 1:ndraw
                if display && mod(i,100) == 0, disp(i); end
                beta = sampling_factor_loading_normal(obj.y, factor, obj.beta_prior, obj.V_prior, diag(sigma_squared));   % sampling static factor loading
                sigma_squared = sampling_sigma_squared(obj.y-factor*beta', obj.gamma_prior, obj.delta_prior);  % sampling static sigma_squared
                for j = 1:obj.k
                    theta(j) = sampling_beta(factor(1:obj.T-1,j),factor(2:obj.T,j),obj.theta_prior,obj.theta_V_prior,1,'stationary',true);  % sampling factor AR(1) coefficients
                end
                factor = sampling_factor_dynamic(obj.y, beta, theta, diag(sigma_squared), eye(obj.k));  % sampling factors from normal distribution
                if i > burnin
                    obj.sampled_beta(:, :, i-burnin) = beta;
                    obj.sampled_sigma_squared(:, i-burnin) = sigma_squared;
                    obj.sampled_theta(:, i-burnin) = theta;
                    obj.sampled_factor(:, :, i-burnin) = factor;
                end
            end
            if display, disp("Done"); end
        end
    end
end

