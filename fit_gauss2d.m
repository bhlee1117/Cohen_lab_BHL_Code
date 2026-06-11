function [pfit, zhat, R2, gof] = fit_gauss2d(Z, X, Y, useWeights, bounds)
% fit_gauss2d  Fit an axis-aligned 2D Gaussian to an image.
%
% Model:
%   Z(x,y) = A * exp( -((x-x0)^2/(2*sx^2) + (y-y0)^2/(2*sy^2)) ) + z0
%
% Important:
%   x0 and y0 are fixed at the maximum finite pixel of Z.
%
% Usage:
%   pfit = fit_gauss2d(Z)
%   pfit = fit_gauss2d(Z, X, Y)
%   pfit = fit_gauss2d(Z, X, Y, useWeights)
%   pfit = fit_gauss2d(Z, X, Y, useWeights, bounds)
%   [pfit, zhat, R2, gof] = fit_gauss2d(...)
%
% Inputs:
%   Z          - [H x W] numeric image. NaNs are ignored.
%   X          - optional [1 x W] x-axis vector. Default: 1:W.
%   Y          - optional [H x 1] y-axis vector. Default: (1:H)'.
%   useWeights - optional true/false. Default: true.
%   bounds     - optional struct to override lower/upper bounds for sx, sy, z0.
%                Any omitted field uses the automatic default.
%                Fields (all optional):
%                  .sx  [lb, ub]  e.g. [2, 50]
%                  .sy  [lb, ub]  e.g. [2, 50]
%                  .z0  [lb, ub]  e.g. [-100, 100]
%                Example:
%                  bounds.sx = [3, 30];
%                  bounds.sy = [3, 30];
%                  bounds.z0 = [0, Inf];
%
% Output:
%   pfit - [1 x 6] fitted parameter vector:
%            pfit(1) = A    amplitude
%            pfit(2) = x0   fixed x center, maximum-pixel x
%            pfit(3) = y0   fixed y center, maximum-pixel y
%            pfit(4) = sx   x spread, sigma
%            pfit(5) = sy   y spread, sigma
%            pfit(6) = z0   background offset
%
%   zhat - fitted Gaussian image
%   R2   - unweighted R-squared on valid pixels
%   gof  - diagnostics structure

%% ---- 0. Input handling ----

Z = double(Z);
[img_h, img_w] = size(Z);

if nargin < 2 || isempty(X)
    X = 1:img_w;
end

if nargin < 3 || isempty(Y)
    Y = (1:img_h)';
end

if nargin < 4 || isempty(useWeights)
    useWeights = true;
end

if nargin < 5 || isempty(bounds)
    bounds = struct();
end

X = X(:)';    % row vector [1 x W]
Y = Y(:);     % column vector [H x 1]

if numel(X) ~= img_w
    error('fit_gauss2d: numel(X) must match size(Z,2) = %d.', img_w);
end

if numel(Y) ~= img_h
    error('fit_gauss2d: numel(Y) must match size(Z,1) = %d.', img_h);
end

%% ---- 1. Build coordinate matrices and vectorize ----

[Xmat, Ymat] = meshgrid(X, Y);

valid = isfinite(Z) & isfinite(Xmat) & isfinite(Ymat);

if nnz(valid) < 4
    error('fit_gauss2d: Not enough finite pixels to fit Gaussian.');
end

xv = Xmat(valid);
yv = Ymat(valid);
zv = Z(valid);

xy_valid = [xv, yv];

%% ---- 2. Fix center at maximum finite pixel ----

[zmax, max_idx] = max(zv);

x0_fixed = xv(max_idx);
y0_fixed = yv(max_idx);

%% ---- 3. Initial parameter guess ----

z0_init = median(zv);

A_init = zmax - z0_init;

% Since center is fixed at maximum pixel, force bright Gaussian amplitude.
if ~isfinite(A_init) || A_init <= 0
    A_init = max(std(zv), eps);
end

% Use positive signal above baseline to estimate initial spread.
weights_init = max(zv - z0_init, 0);

if sum(weights_init) > 0
    sx_init = sqrt(sum(weights_init .* (xv - x0_fixed).^2) / sum(weights_init));
    sy_init = sqrt(sum(weights_init .* (yv - y0_fixed).^2) / sum(weights_init));
else
    sx_init = (max(X) - min(X)) / 6;
    sy_init = (max(Y) - min(Y)) / 6;
end

%% ---- 4. Coordinate spacing and bounds ----

dx_vals = diff(unique(X));
dy_vals = diff(unique(Y));

if isempty(dx_vals)
    dx = 1;
else
    dx = median(abs(dx_vals));
end

if isempty(dy_vals)
    dy = 1;
else
    dy = median(abs(dy_vals));
end

if ~isfinite(dx) || dx <= 0
    dx = 1;
end

if ~isfinite(dy) || dy <= 0
    dy = 1;
end

x_range = max(X) - min(X);
y_range = max(Y) - min(Y);

if x_range <= 0
    x_range = dx;
end

if y_range <= 0
    y_range = dy;
end

sx_min = dx / 2;
sy_min = dy / 2;

sx_max = 2 * x_range;
sy_max = 2 * y_range;

sx_init = min(max(sx_init, sx_min), sx_max);
sy_init = min(max(sy_init, sy_min), sy_max);

%% ---- 5. Define model with fixed center ----
% q = [A, sx, sy, z0]
%
% pfit returned later is:
% pfit = [A, x0_fixed, y0_fixed, sx, sy, z0]

gauss2d_fixed_center = @(q, xy) ...
    q(1) .* exp( ...
    -( ...
    (xy(:,1) - x0_fixed).^2 ./ (2 * q(2)^2) + ...
    (xy(:,2) - y0_fixed).^2 ./ (2 * q(3)^2) ...
    ) ...
    ) + q(4);

q0 = [A_init, sx_init, sy_init, z0_init];

% Bright Gaussian peak because center is fixed at maximum pixel.
% Apply user-supplied bounds where provided, else use auto defaults.
if isfield(bounds, 'sx') && numel(bounds.sx) == 2
    sx_lb = bounds.sx(1);   sx_ub = bounds.sx(2);
else
    sx_lb = sx_min;         sx_ub = sx_max;
end

if isfield(bounds, 'sy') && numel(bounds.sy) == 2
    sy_lb = bounds.sy(1);   sy_ub = bounds.sy(2);
else
    sy_lb = sy_min;         sy_ub = sy_max;
end

if isfield(bounds, 'z0') && numel(bounds.z0) == 2
    z0_lb = bounds.z0(1);   z0_ub = bounds.z0(2);
else
    z0_lb = -Inf;           z0_ub = Inf;
end

% Clamp initial guesses to respect user bounds
sx_init = min(max(sx_init, sx_lb), sx_ub);
sy_init = min(max(sy_init, sy_lb), sy_ub);
z0_init = min(max(z0_init, z0_lb), z0_ub);

lb = [0,   sx_lb, sy_lb, z0_lb];
ub = [Inf, sx_ub, sy_ub, z0_ub];

%% ---- 6. Optional weighted least-squares ----

if useWeights
    % Intensity-based weights.
    % Pixels above baseline get more weight, but all pixels still contribute.
    w = max(zv - z0_init, 0);

    if max(w) > 0
        w = w / max(w);

        % Weight floor prevents the background/tails from being completely ignored.
        weightFloor = 0.05;
        w = weightFloor + (1 - weightFloor) * w;
    else
        w = ones(size(zv));
    end
else
    w = ones(size(zv));
end

sqrtw = sqrt(w);

weighted_model = @(q, xy) sqrtw .* gauss2d_fixed_center(q, xy);
weighted_data  = sqrtw .* zv;

%% ---- 7. Fit ----

if exist('lsqcurvefit', 'file') ~= 2
    error('fit_gauss2d requires lsqcurvefit from the Optimization Toolbox.');
end

opts = optimoptions('lsqcurvefit', ...
    'Display', 'off', ...
    'MaxFunctionEvaluations', 1e5, ...
    'MaxIterations', 1e4);

[qfit, resnorm, residual_weighted, exitflag, output] = ...
    lsqcurvefit(weighted_model, q0, xy_valid, weighted_data, lb, ub, opts);

%% ---- 8. Convert qfit to pfit ----

A_fit  = qfit(1);
sx_fit = qfit(2);
sy_fit = qfit(3);
z0_fit = qfit(4);

pfit = [A_fit, x0_fixed, y0_fixed, sx_fit, sy_fit, z0_fit];

%% ---- 9. Evaluate fitted image ----

zhat_vec = gauss2d_fixed_center(qfit, [Xmat(:), Ymat(:)]);
zhat = reshape(zhat_vec, img_h, img_w);

resid = Z - zhat;

zfit_valid = gauss2d_fixed_center(qfit, xy_valid);
resid_valid = zv - zfit_valid;

SSR = sum(resid_valid.^2);
SST = sum((zv - mean(zv)).^2);

if SST > eps
    R2 = 1 - SSR / SST;
else
    R2 = NaN;
end

%% ---- 10. Weighted R2, if desired ----

zbar_w = sum(w .* zv) / sum(w);
SSR_w = sum(w .* resid_valid.^2);
SST_w = sum(w .* (zv - zbar_w).^2);

if SST_w > eps
    R2_weighted = 1 - SSR_w / SST_w;
else
    R2_weighted = NaN;
end

%% ---- 11. Package diagnostics ----

gof = struct();

gof.FWHM_x = 2 * sqrt(2 * log(2)) * sx_fit;
gof.FWHM_y = 2 * sqrt(2 * log(2)) * sy_fit;

gof.residual = resid;

gof.R2 = R2;
gof.R2_weighted = R2_weighted;

gof.resnorm = resnorm;
gof.residual_weighted = residual_weighted;
gof.exitflag = exitflag;
gof.output = output;

gof.center_fixed = true;
gof.x0_fixed = x0_fixed;
gof.y0_fixed = y0_fixed;
gof.max_pixel_value = zmax;

gof.useWeights = useWeights;

weights_image = NaN(size(Z));
weights_image(valid) = w;
gof.weights = weights_image;

if exitflag <= 0
    warning('fit_gauss2d: lsqcurvefit may not have converged. exitflag = %d.', exitflag);
end

end