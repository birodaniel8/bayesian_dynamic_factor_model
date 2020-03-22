%% Generate data:
x = bfmdgp();
%% Run normal linear model estimation:
k = 3;
priors = struct([]);
model = BDFM(x, k, priors);
n_draw = 1500;
n_burnin = 1000;
[sampled_beta, sampled_sigma_squared, sampled_factor] = model.estimate_normal_model(n_draw, n_burnin);
%% Estimated values:
mean(sampled_beta,3)
mean(sampled_sigma_squared,2)
%% Plot estimated factors:
plot(mean(sampled_factor,3))
%% Plot beta(1,1) histogram:
hist(squeeze(sampled_beta(1,1,:)))
%% Significance of the variables:
alpha = 0.95;
sign(quantile(sampled_beta,(1-alpha)/2,3)) .* sign(quantile(sampled_beta,1-(1-alpha)/2,3))


%% Import data:
load('currency_data.mat')
%% Run normal linear model estimation:
k = 3;
priors = struct([]);
model = BDFM(x, k, priors);
n_draw = 1500;
n_burnin = 1000;
[sampled_beta, sampled_sigma_squared, sampled_factor] = model.estimate_normal_model(n_draw, n_burnin);
%% Estimated values:
mean(sampled_beta,3)
mean(sampled_sigma_squared,2)
%% Plot estimated factors:
plot(mean(sampled_factor,3))
%% Plot beta(1,1) histogram:
hist(squeeze(sampled_beta(1,1,:)))
%% Significance of the variables:
alpha = 0.99;
sign(quantile(sampled_beta,(1-alpha)/2,3)) .* sign(quantile(sampled_beta,1-(1-alpha)/2,3))