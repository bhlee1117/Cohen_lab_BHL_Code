function plot_StackLine(dataMatrix, xAxis, layerSpacing, fillAlpha)
    % Function to create a ridgeline plot from a data matrix
    % Inputs:
    %   dataMatrix   - n x T matrix (data to be plotted)
    %   xAxis        - 1 x T vector (X-axis values)
    %   layerSpacing - Scalar, spacing between layers (default: 1)
    %   lineAlpha    - Scalar, transparency of the lines (default: 0.8)
    %   fillAlpha    - Scalar, transparency of the fill (default: 0.5)

    % Default values
    if nargin < 3, layerSpacing = 1; end
    if nargin < 4, fillAlpha = 1; end

    % Validate dimensions
    [n, T] = size(dataMatrix);
    if numel(xAxis) ~= T
        error('xAxis must have the same number of elements as columns in dataMatrix.');
    end

    % Set up the figure
    clf;
    hold on;

    % Define colormap
    cmap = turbo(n); % Color map with `n` layers

    % Loop through each layer and plot
    for i = [n:-1:1]
        % Offset the layer by its position
        yOffset = (i - 1) * layerSpacing;

        % Remove NaN values
        validIndices = ~isnan(dataMatrix(i, :));
        xValid = xAxis(validIndices);
        zValid = dataMatrix(i, validIndices);

        % Fill the area under the curve
        fill([xValid, fliplr(xValid)], ...
             [zValid + yOffset, yOffset * ones(size(zValid))], ...
             cmap(i, :), 'FaceAlpha', fillAlpha, 'EdgeColor', 'none');

        % Add the line on top
        plot(xValid, zValid + yOffset, 'Color', [0 0 0], 'LineWidth', 1.2);
    end

    % Adjust axes
    axis tight;
    grid on;
    hold off;
end