function BreakYAxis(X, Y, start, stop, width, colormap)
% BreakYAxis creates a folded Y-axis for multiple traces with custom colors.
% 
% Inputs:
%   X        - 1xT vector of X-axis values (time or other independent variable)
%   Y        - NxT matrix of Y-axis values (N traces, T time points)
%   start    - Lower bound of the break on the Y-axis
%   stop     - Upper bound of the break on the Y-axis
%   width    - Gap width (in Y-axis units) for the visual break
%   colormap - Nx3 matrix specifying RGB colors for each trace
%              (or use a MATLAB colormap function like `jet(N)` for N traces)
%
% Example:
%   T = 100; N = 5;
%   X = linspace(0, 10, T);
%   Y = rand(N, T); % Random data for 5 traces
%   colors = jet(N); % Use jet colormap
%   BreakYAxis(X, Y, 0.2, 0.8, 0.05, colors);

% Validate inputs
if start >= stop
    error('The "start" value must be less than the "stop" value.');
end
if size(X, 2) ~= size(Y, 2)
    error('The size of X must match the number of columns in Y.');
end
if size(Y, 1) ~= size(colormap, 1)
    error('The number of rows in Y must match the number of colors in the colormap.');
end

% Remove data in the break range
YModified = Y;
YModified(Y > start & Y < stop) = NaN;

% Shift data above the break downward
YShifted = YModified;
YShifted(Y > stop) = YShifted(Y > stop) - (stop - start - width);

% Plot the data with custom colors
hold on;
for i = 1:size(YShifted, 1)
    plot(X, YShifted(i, :), 'LineWidth', 1.5, 'Color', colormap(i, :));
end

% Add the folded axis indicators ("//")
xLimits = get(gca, 'XLim');
xSpan = diff(xLimits) * 0.02; % Length of the slashes relative to X-axis range
plot([xLimits(1), xLimits(1) + xSpan], [start, start], 'k-', 'LineWidth', 1.5);
plot([xLimits(2) - xSpan, xLimits(2)], [start, start], 'k-', 'LineWidth', 1.5);
plot([xLimits(1), xLimits(1) + xSpan], [start - (stop - start - width), start - (stop - start - width)], 'k-', 'LineWidth', 1.5);
plot([xLimits(2) - xSpan, xLimits(2)], [start - (stop - start - width), start - (stop - start - width)], 'k-', 'LineWidth', 1.5);

% Adjust the Y-axis limits
ylim([min(YShifted(:)) - 0.1, max(YShifted(:)) + 0.1]);

% Adjust the Y-ticks to reflect the folded axis
yTicks = get(gca, 'YTick');
newTicks = yTicks(yTicks <= start | yTicks >= stop) - (yTicks > stop) * (stop - start - width);
set(gca, 'YTick', newTicks);

% Customize the Y-tick labels
yTickLabels = arrayfun(@(yt) sprintf('%.2f', yt), newTicks, 'UniformOutput', false);
set(gca, 'YTickLabel', yTickLabels);

xlabel('X-axis');
ylabel('Y-axis');
title('Folded Y-axis with Break');
hold off;
end