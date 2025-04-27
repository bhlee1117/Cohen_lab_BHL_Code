function [binned_data, x_edges, y_edges] = binning_data2D(X, Y, values, x_edges, y_edges)
    % BIN2D Performs 2D binning of data using custom bin edges
    %   [binned_data, x_edges, y_edges] = bin2d(X, Y, values, x_edges, y_edges)
    %
    %   X, Y       - Data coordinates (vectors of the same size)
    %   values     - Values to be binned (same size as X and Y)
    %   x_edges    - Custom bin edges along X-axis (vector)
    %   y_edges    - Custom bin edges along Y-axis (vector)
    %
    %   binned_data - 2D matrix of binned values (sum of values in each bin)
    %   x_edges, y_edges - Bin edges for X and Y (same as input)

    if numel(X) ~= numel(Y) || numel(X) ~= numel(values)
        error('X, Y, and values must have the same length.');
    end

    % Assign data to bins using discretization
    x_idx = discretize(X, x_edges);
    y_idx = discretize(Y, y_edges);

    % Remove out-of-bounds indices
    valid = ~isnan(x_idx) & ~isnan(y_idx);
    x_idx = x_idx(valid);
    y_idx = y_idx(valid);
    values = values(valid);

    % Initialize binned data matrix
    binned_data = zeros(numel(y_edges) - 1, numel(x_edges) - 1);

    % Accumulate values into bins
    for i = 1:numel(values)
        binned_data(y_idx(i), x_idx(i)) = binned_data(y_idx(i), x_idx(i)) + values(i);
    end
end
