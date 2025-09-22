function h = errorbar_shade(x, y, error, cmap)
% errorbar_shade - Plots a line with a shaded error region, robust to NaNs
%
% Syntax:
%   h = errorbar_shade(x, y, error)
%   h = errorbar_shade(x, y, error, cmap)
%
% Inputs:
%   x     - x values (vector)
%   y     - y values (vector)
%   error - error margins (vector, same size as y), can contain NaNs
%   cmap  - color (optional, RGB triplet)
%
% Output:
%   h - handle to the main plot line

    % Ensure row vectors
    x = x(:)';
    y = y(:)';
    error = error(:)';

    % Validate input sizes
    if ~isequal(size(x), size(y)) || ~isequal(size(y), size(error))
        error('x, y, and error must be vectors of the same size.');
    end

    % Default color
    if nargin < 4 || isempty(cmap)
        cmap = distinguishable_colors(1);
        cmap = cmap(1, :);
    end

    % Create upper/lower curves, ignoring NaNs
    upper = y + error;
    lower = y - error;

    % Find finite (non-NaN) segments to shade
    valid = isfinite(upper) & isfinite(lower) & isfinite(x);
    if any(valid)
        fill([x(valid), fliplr(x(valid))], ...
             [upper(valid), fliplr(lower(valid))], ...
             cmap, 'FaceAlpha', 0.3, 'LineStyle', 'none');
        hold on;
    end

    % Plot central line (will handle its own NaNs)
    h = plot(x, y, 'Color', cmap, 'LineWidth', 2);
end
