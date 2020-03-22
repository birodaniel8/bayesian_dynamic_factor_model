%% Generate data:
x = bfmdgp();
%% Run normal linear model estimation:
k = 3;
priors = struct([]);
model = BDFM_nlm(x, k, priors);
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


%% Import data:
load('currency_data.mat')
%% Run normal linear model estimation:
k = 2;
priors = struct([]);
model = BDFM_nlm(x, k, priors);
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