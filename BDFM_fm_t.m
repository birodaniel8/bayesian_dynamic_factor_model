classdef BDFM_fm_t < BDFM
    %BDFM_fm_t Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        beta_prior, V_prior, gamma_prior, delta_prior, v_prior
        mh_sigma_squared
        sampled_beta, sampled_sigma_squared, sampled_df, sampled_scale_parameter, sampled_factor
    end
    
    methods
        function obj = BDFM_fm_t(y, k, priors)
            %BDFM_nlm Construct an instance of this class
            obj@BDFM(y, k);
            if isfield(priors,'beta'), obj.beta_prior = priors.beta; else obj.beta_prior = 0; end  % prior mean of the beta coefficients in the normal model
            if isfield(priors,'V'), obj.V_prior = priors.V; else obj.V_prior = 1; end  % prior variance of the beta coefficients in the normal model
            if isfield(priors,'gamma'), obj.gamma_prior = priors.gamma; else obj.gamma_prior = 3/2; end  % prior confidence in the variance in the normal model
            if isfield(priors,'delta'), obj.delta_prior = priors.delta; else obj.delta_prior = 0.01/2; end  % prior mean of the variance in the normal model
            if isfield(priors,'v'), obj.v_prior = priors.v; else obj.v_prior = 30; end  % prior degfree of freedom parameter
            if isfield(priors,'mh_sigma_squared'), obj.mh_sigma_squared = priors.mh_sigma_squared; else obj.mh_sigma_squared = 1; end  % random walk Metropolis-Hastings algorithm variance parameters
        end
        
        function obj = estimate(obj, ndraw, burnin, display)
            if nargin == 3, display = true; end
            % create containers:
            obj.sampled_beta = zeros(obj.m, obj.k, ndraw-burnin);
            obj.sampled_sigma_squared = zeros(obj.m, ndraw-burnin);
            obj.sampled_df = zeros(obj.m, ndraw-burnin);
            obj.sampled_scale_parameter = zeros(obj.T, obj.m, ndraw-burnin);
            obj.sampled_factor = zeros(obj.T, obj.k, ndraw-burnin);
            
            % initial values:
            sigma_squared = repmat(obj.delta_prior, obj.m, 1);
            factor = factor_initialize(obj.y,obj.k);  % revisit this function (this it not my own function)
            df = repmat(obj.v_prior, obj.m, 1);
            scale_parameter = ones(obj.T, obj.m);
            
            % sampling:
            for i = 1:ndraw
                if display && mod(i,100) == 0, disp(i); end
                beta = sampling_factor_loading_normal(obj.y, factor, obj.beta_prior, obj.V_prior, repmat(sigma_squared,1,obj.T) .* (1 ./ scale_parameter)');   % sampling static factor loading
                sigma_squared = sampling_sigma_squared((obj.y-factor*beta') .* sqrt(scale_parameter), obj.gamma_prior, obj.delta_prior);  % sampling static sigma_squared
                df = sampling_degree_of_freedom(scale_parameter,df,obj.v_prior,obj.mh_sigma_squared);  % sampling the degree of freedom parameter
                scale_parameter = sampling_mixture_scale_parameter(obj.y-factor*beta',sigma_squared,df);  % sampling the scale mixture parameter
                factor = sampling_factor_normal(obj.y, beta, repmat(sigma_squared,1,obj.T) .* (1 ./ scale_parameter)');  % sampling factors from normal distribution
                if i > burnin
                    obj.sampled_beta(:, :, i-burnin) = beta;
                    obj.sampled_sigma_squared(:, i-burnin) = sigma_squared;
                    obj.sampled_df(:, i-burnin) = df;
                    obj.sampled_scale_parameter(:, :, i-burnin) = scale_parameter;
                    obj.sampled_factor(:, :, i-burnin) = factor;
                end
            end
            if display, disp("Done"); end
        end
    end
end

