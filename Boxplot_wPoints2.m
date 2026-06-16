function p_values = Boxplot_wPoints2(X, faceColor)
% p_values = Boxplot_wPoints2(X, faceColor)
% Boxplot_wPoints2 - Plots box plots with beeswarm or jittered points.
% Supports both unpaired (cell array) and paired (numeric matrix) data.
%
%   - Cell array input → unpaired, ranksum.
%   - Numeric matrix input → paired, ttest.
%
%   Also connects points with lines for paired input.
markersize=10;
isPaired = isnumeric(X) && ismatrix(X);
if isPaired
    numGroups = size(X,2);
    data = X;
elseif iscell(X)
    numGroups = numel(X);
    data = X;
else
    error('Input must be either a cell array (unpaired) or numeric matrix (paired).');
end

% --- Handle faceColor ---
if ischar(faceColor) || isstring(faceColor)
    colorMap = containers.Map({'r','g','b','c','m','y','k','w'}, ...
        {[1 0 0], [0 1 0], [0 0 1], [0 1 1], [1 0 1], [1 1 0], [0 0 0], [1 1 1]});
    if isKey(colorMap, lower(faceColor))
        faceColor = repmat(colorMap(lower(faceColor)), numGroups, 1);
    else
        error('Unknown color specification.');
    end
end
if size(faceColor,1) == 1
    faceColor = repmat(faceColor, numGroups, 1);
elseif size(faceColor,1) ~= numGroups
    error('faceColor must match number of groups.');
end

hold on;

% --- Plot boxplots ---
for i = 1:numGroups
    if isPaired
        values = data(:,i);
    else
        values = data{i};
    end
    values = values(~isnan(values));
    
    if exist('boxchart', 'file')
        boxchart(i * ones(size(values)), values, ...
            'BoxFaceColor', faceColor(i,:), 'MarkerStyle', 'none');
    else
        boxplot(values, 'Positions', i, 'Colors', 'k', 'Symbol', '');
    end
end

% --- Beeswarm / jittered points ---
for i = 1:numGroups
    if isPaired
        values = data(:,i);
    else
        values = data{i};
    end
    xVals = repmat(i, size(values));
    
    if exist('swarmchart', 'file')
        swarmchart(xVals, values, markersize, faceColor(i,:), 'filled', ...
            'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', 'k', ...
            'XJitter', 'randn', 'XJitterWidth', 0.1);
    else
        x_jittered = i + 0.1 * (rand(size(values)) - 0.5);
        scatter(x_jittered, values, markersize, faceColor(i,:), 'filled', ...
            'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', 'k');
    end
end

% --- Connect paired points ---
if isPaired
    for i = 1:size(data,1)
        validCols = find(~isnan(data(i,:)));
        if numel(validCols) >= 2
            plot(validCols, data(i,validCols), '-k', 'LineWidth', 0.5, 'Color', [0.5 0.5 0.5 0.5]);
        end
    end
end

% --- Statistical tests ---
p_values = nan(numGroups);
for i = 1:numGroups
    for j = i+1:numGroups
        if isPaired
            pairedData = [data(:,i), data(:,j)];
            mask = all(~isnan(pairedData),2);
            [~, p] = ttest(pairedData(mask,1), pairedData(mask,2));
        else
            xi = data{i}(~isnan(data{i}));
            xj = data{j}(~isnan(data{j}));
            if ~isempty(xi) & ~isempty(xj)
            p = ranksum(xi, xj);
            else
            p = NaN;
            end
        end
        p_values(i,j) = p;
        p_values(j,i) = p;
    end
end

% --- Axis formatting ---
set(gca, 'XTick', 1:numGroups, ...
    'XTickLabel', arrayfun(@(i) sprintf('Group %d', i), 1:numGroups, 'UniformOutput', false));
ylabel('Value');
box on;
hold off;
end