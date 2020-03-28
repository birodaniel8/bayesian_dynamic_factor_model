%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Estimate FM on simulated data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate data:
x = bfmdgp();
%% Run normal linear model estimation:
k = 3;
priors = struct();
model = BDFM_fm(x, k, priors);
n_draw = 1500;
n_burnin = 1000;
model.estimate(n_draw, n_burnin);
%% Estimated values:
mean(model.sampled_beta,3)
mean(model.sampled_sigma_squared,2)'
%% Plot estimated factors:
plot(mean(model.sampled_factor,3))
%% Plot beta(1,1) histogram:
hist(squeeze(model.sampled_beta(1,1,:)))
%% Significance of the variables:
alpha = 0.95;
sign(quantile(model.sampled_beta,(1-alpha)/2,3)) .* sign(quantile(model.sampled_beta,1-(1-alpha)/2,3))
%% estimates with t distribution (df=30 as prior)
model_t = BDFM_fm_t(x, k, priors);
model_t.estimate(n_draw, n_burnin);
%% Estimated values:
mean(model_t.sampled_beta,3)
mean(model_t.sampled_sigma_squared,2)'
%% Significance of the variables:
alpha = 0.95;
sign(quantile(model_t.sampled_beta,(1-alpha)/2,3)) .* sign(quantile(model_t.sampled_beta,1-(1-alpha)/2,3))

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Estimate FM on simulated data with t distribution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate data:
x = bfmdgp_t();
%% Run normal linear model estimation:
k = 3;
priors = struct();
model = BDFM_fm(x, k, priors);
n_draw = 1500;
n_burnin = 1000;
model.estimate(n_draw, n_burnin);
%% Estimated values:
mean(model.sampled_beta,3)
mean(model.sampled_sigma_squared,2)'
%% Significance of the variables:
alpha = 0.95;
sign(quantile(model.sampled_beta,(1-alpha)/2,3)) .* sign(quantile(model.sampled_beta,1-(1-alpha)/2,3))
%% estimates with t distribution (df=5 as prior)
priors.v_prior = 5;
model_t = BDFM_fm_t(x, k, priors);
model_t.estimate(n_draw, n_burnin);
%% Estimated values:
mean(model_t.sampled_beta,3)
mean(model_t.sampled_sigma_squared,2)'
%% Significance of the variables:
alpha = 0.95;
sign(quantile(model_t.sampled_beta,(1-alpha)/2,3)) .* sign(quantile(model_t.sampled_beta,1-(1-alpha)/2,3))




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Estimate DFM on simulated data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate data:
x = bfmdgp();
%% Run normal linear model estimation:
k = 3;
priors = struct();
model = BDFM_dfm(x, k, priors);
n_draw = 1500;
n_burnin = 1000;
model.estimate(n_draw, n_burnin);
%% Estimated values:
mean(model.sampled_beta,3)
mean(model.sampled_sigma_squared,2)'
mean(model.sampled_theta,2)'
%% Plot estimated factors:
plot(mean(model.sampled_factor,3))
%% Plot beta(1,1) histogram:
hist(squeeze(model.sampled_beta(1,1,:)))
%% Significance of the variables:
alpha = 0.95;
sign(quantile(model.sampled_beta,(1-alpha)/2,3)) .* sign(quantile(model.sampled_beta,1-(1-alpha)/2,3))
sign(quantile(model.sampled_theta,(1-alpha)/2,2)) .* sign(quantile(model.sampled_theta,1-(1-alpha)/2,2))

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Estimate FM on currency data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Import data:
load('currency_data.mat')
%% Run normal linear model estimation:
k = 2;
priors = struct();
model = BDFM_fm(x, k, priors);
n_draw = 1500;
n_burnin = 1000;
model.estimate(n_draw, n_burnin);
%% Estimated values:
mean(model.sampled_beta,3)
mean(model.sampled_sigma_squared,2)'
%% Plot estimated factors:
plot(mean(model.sampled_factor,3))
%% Plot beta(1,1) histogram:
hist(squeeze(model.sampled_beta(1,1,:)))
%% Significance of the variables:
alpha = 0.95;
sign(quantile(model.sampled_beta,(1-alpha)/2,3)) .* sign(quantile(model.sampled_beta,1-(1-alpha)/2,3))