function [fit_curves, cFitted, R2] = expfit_wBd(x, y, xv, initialGuess, lowerBounds, upperBounds, weights)
% [fit_curves, cFitted, R2] = expfit_wBd(x, y, xv, initialGuess, lowerBounds, upperBounds)
% [fit_curves, cFitted, R2] = expfit_wBd(x, y, xv, initialGuess, lowerBounds, upperBounds, weights)
% 2024.10.28 Byung Hun Lee, exponential fitting
% weights (optional): vector of non-negative weights, same length as x
% Define the exponential model function: y = a*exp(b*x)
modelFun = @(c, x) c(1) * exp(x/c(2));

if nargin < 7 || isempty(weights)
    weights = ones(size(x));
end
weights = weights(:) / sum(weights);

% Perform fitting with lsqcurvefit using weighted residuals
weightedResidFun = @(c) sqrt(weights) .* (modelFun(c, x(:)) - y(:));
options = optimoptions('lsqnonlin', 'Display', 'off');
cFitted = lsqnonlin(weightedResidFun, initialGuess, lowerBounds, upperBounds, options);

fit_curves = modelFun(cFitted, xv);

% Weighted R^2 computed on the original data points (x, y)
y_pred   = modelFun(cFitted, x(:));
SS_res   = sum(weights .* (y(:) - y_pred).^2);
SS_tot   = sum(weights .* (y(:) - sum(weights .* y(:))).^2);
R2       = 1 - SS_res / SS_tot;
end