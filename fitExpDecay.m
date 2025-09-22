function [y_fitc, params, R2] = fitExpDecay(x, y, xc)
% fitExpDecay: Fits data to a * exp(-x / b) with heavier weight on first 3 points
% [y_fitc, params, R2] = fitExpDecay(x, y, xc)
%
% Inputs:
%   x - Independent variable (vector)
%   y - Dependent variable (vector)
%   xc - Points to evaluate the fitted curve
%
% Outputs:
%   y_fitc - Fitted values at xc
%   params - [a, b] coefficients
%   R2 - Coefficient of determination (on original x)

    % Remove NaNs
    valid = ~isnan(x) & ~isnan(y);
    x = x(valid);
    y = y(valid);

    % Initial guess
    a0 = max(y) - min(y);
    b0 = (max(x) - min(x)) / 2;
    p0 = [a0, b0];

    % Weights: boost first 3 points
    weights = ones(size(y));
    weights(1:min(3,end)) = 5;  % give 5x weight to first 3 points

    % Weighted residual
    model = @(p, x) p(1) * exp(-x / p(2));
    residual = @(p) sqrt(weights(:)) .* (model(p, x) - y(:));

    % Optimization
    lb = [-Inf, 30];
    ub = [Inf, 2000];
    opts = optimoptions('lsqnonlin', 'Display', 'off');
    params = lsqnonlin(residual, p0, lb, ub, opts);

    % Output
    y_fitc = model(params, xc);
    y_fit = model(params, x);

    % R² (on original data)
    SSres = sum((y - y_fit).^2);
    SStot = sum((y - mean(y)).^2);
    R2 = 1 - SSres / SStot;
end
