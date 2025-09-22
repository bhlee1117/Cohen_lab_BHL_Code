function [p, r_squared, residual, kept_indices, polyfunc] = fit_poly(x, y, w, max_iterations, threshold_factor)
% fit_poly Performs iterative weighted linear polynomial fitting
%   [p, r_squared, residual, kept_indices] = fit_poly(x, y, w, max_iterations, threshold_factor)
%   fits a linear polynomial to data points (x,y) with weights w, iteratively
%   removing outliers based on residual threshold.
%
%   Inputs:
%       x - vector of x-coordinates
%       y - vector of y-coordinates
%       w - vector of weights (same length as x and y)
%       max_iterations - maximum number of iterations for outlier removal
%       threshold_factor - multiplier for standard deviation to define outliers
%                         (e.g., 2 means residuals > 2*std(residual) are removed)
%
%   Outputs:
%       p - coefficients of the linear polynomial [slope, intercept]
%       r_squared - coefficient of determination
%       residual - residual vector for kept points
%       kept_indices - indices of data points retained after outlier removal

%   Inputs:
%       x - vector of x-coordinates
%       y - vector of y-coordinates
%       w - vector of weights (same length as x and y)
%       max_iterations - maximum number of iterations for outlier removal
%       threshold_factor - multiplier for standard deviation to define outliers
%                         (e.g., 2 means residuals > 2*std(residual) are removed)
%
%   Outputs:
%       p - coefficients of the linear polynomial [slope, intercept]
%       r_squared - coefficient of determination
%       residual - residual vector for kept points (NaN for removed points)
%       kept_indices - indices of data points retained after outlier removal

% Input validation
if nargin < 3
    w = ones(size(x)); % Default to equal weights if none provided
end
if nargin < 4
    max_iterations = 10; % Default maximum iterations
end
if nargin < 5
    threshold_factor = 2; % Default threshold (2 standard deviations)
end
if length(x) ~= length(y) || length(x) ~= length(w)
    error('Input vectors x, y, and w must have the same length');
end
if any(w < 0)
    error('Weights must be non-negative');
end
if max_iterations < 1
    error('max_iterations must be at least 1');
end

% Ensure column vectors
x = x(:);
y = y(:);
w = w(:);

% Initialize
kept_indices = true(size(x));
current_x = x;
current_y = y;
current_w = w;

% Create fit options with weights
fo = fitoptions('Method', 'LinearLeastSquares', 'Weights', current_w, 'Lower', [-Inf, -Inf], 'Upper', [Inf, Inf]);

for iter = 1:max_iterations
    % Update fit options with current weights
    fo.Weights = current_w;
    
    % Perform weighted linear regression using fit
    f = fit(current_x, current_y, 'poly1', fo);
    
    % Extract coefficients (slope, intercept)
    p = [f.p1, f.p2];
    
    % Compute fitted values and residuals
    y_fit = f(current_x);
    residual = current_y - y_fit;
    
    % Compute standard deviation of residuals
    residual_std = std(residual);
    
    % Identify outliers (residuals exceeding threshold_factor * std)
    keep_points = abs(residual) <= threshold_factor * residual_std;
    
    % If no points are removed or too few points remain, stop
    if all(keep_points) || sum(keep_points) < 2
        break;
    end
    
    % Update data by removing outliers
    current_x = current_x(keep_points);
    current_y = current_y(keep_points);
    current_w = current_w(keep_points);
    
    % Update kept_indices
    temp_indices = find(kept_indices);
    kept_indices(temp_indices(~keep_points)) = false;
end

% Final fit with remaining points
fo.Weights = current_w;
f = fit(current_x, current_y, 'poly1', fo);
p = [f.p1, f.p2];
y_fit = f(current_x);
residual = current_y - y_fit;

% Compute R-squared
y_mean = sum(current_w.*current_y)/sum(current_w);
SS_tot = sum(current_w.*(current_y - y_mean).^2);
SS_res = sum(current_w.*residual.^2);
r_squared = 1 - SS_res/SS_tot;

% Return residuals for all original points, setting removed points to NaN
full_residual = nan(size(x));
full_residual(kept_indices) = residual;

residual = full_residual;
kept_indices = find(kept_indices); % Convert to indices

polyfunc=@(x) p(1)*x+p(2);
end