function estimated_model = estimate_BDFM(y, k, priors, model, ndraw, burnin, display)
%estimate_BDFM Summary of this function goes here
%   Detailed explanation goes here

if nargin == 6, display = true; end
if nargin >= 6
    if model == "NLM"
        estimated_model = BDFM_nlm(y, k, priors);
        estimated_model.estimate(ndraw, burnin, display);
    else
        estimated_model = -1;
    end
else
    estimated_model = -1;
end

