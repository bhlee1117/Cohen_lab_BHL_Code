function [y_fitc, params, R2] = fitExpDecayC(x, y, xc)
% fitExpDecay: Fits data to a * exp(-x / b) + c
%
% Inputs:
%   x - Independent variable (e.g., time), vector
%   y - Dependent variable (e.g., signal), vector
%
% Outputs:
%   y_fit - Fitted y values
%   params - [a, b, c] coefficients
%   R2 - Coefficient of determination

    if length(x) ~= length(y)
        error('x and y must be the same length.');
    end

    validind=find(~isnan(x) & ~isnan(y));
    x=x(validind);
    y=y(validind);
    % Initial parameter guess: [a, b, c]
    a0 = max(y) - min(y);
    b0 = (max(x) - min(x)) / 2;
    c0 = min(y);
    p0 = [a0, b0, c0];
    %p0 = [a0, b0];

    % Fit model using nonlinear least squares
    %model = @(p, x) p(1) * exp(-x / p(2));
    model = @(p, x) p(1) * exp(-x / p(2)) + p(3);
    opts = optimoptions('lsqcurvefit', 'Display', 'off');
    lb = [-Inf, 0, -Inf];  % b must be positive
    ub = [Inf, Inf, Inf];
    % lb = [-Inf, 0];  % b must be positive
    % ub = [Inf, 2000];
    % lb = [ 0];  % b must be positive
    % ub = [2000];

    params = lsqcurvefit(model, p0, x, y, lb, ub, opts);

    y_fitc = model(params, xc);
    y_fit = model(params, x);

    % Compute R²
    SSres = sum((y - y_fit).^2);
    SStot = sum((y - mean(y)).^2);
    R2 = 1 - SSres / SStot;
end
