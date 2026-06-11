function p_values = Violin_wPoints(X, faceColor)
% Violin_wPoints  Plots violin plots with overlaid data points and stats
%
%   p_values = Violin_wPoints(X, faceColor) uses your custom Violin class
%   to render violin plots for each group in X and performs pairwise
%   independent t-tests, returning a matrix of p-values.
%
%   Inputs:
%       X         - Either:
%                     (a) Cell array {X{1}, ..., X{N}} of vectors
%                     (b) Numeric matrix (N × R), each column is a group
%       faceColor - Either:
%                     - a single RGB or color char
%                     - an N×3 RGB matrix
%
%   Output:
%       p_values  - Symmetric matrix of p-values from t-tests

% --- Convert matrix to cell array if needed ---
if isnumeric(X)
    X = num2cell(X, 1);
elseif iscell(X)
    if ~all(cellfun(@isnumeric, X))
        error('All elements in cell array X must be numeric vectors.');
    end
else
    error('Input X must be a numeric matrix or a cell array of vectors.');
end

numGroups = numel(X);

% --- Normalize faceColor ---
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
    error('faceColor must be 1x3 or Nx3 RGB matrix.');
end

% --- Plot ---
hold on;
for i = 1:numGroups
    Violin({X{i}}, i, ...
        'ViolinColor', {faceColor(i,:)}, ...
        'EdgeColor', [0.3 0.3 0.3], ...
        'ShowMean', false, ...
        'ShowNotches', false, ...
        'ShowBox', true, ...
        'ShowMedian', true, ...
        'ShowData', true, ...
        'HalfViolin', 'full');
end
xlim([0.5, numGroups + 0.5]);
set(gca, 'XTick', 1:numGroups, ...
         'XTickLabel', arrayfun(@(i) sprintf('Group %d', i), 1:numGroups, 'UniformOutput', false));
ylabel('Value');
hold off;

% --- Statistics ---
p_values = nan(numGroups);
for i = 1:numGroups
    for j = i+1:numGroups
        xi = X{i}(~isnan(X{i}));
        xj = X{j}(~isnan(X{j}));
        if numel(xi) > 1 && numel(xj) > 1
            %[~, p] = ttest2(xi, xj);
            [p] = ranksum(xi, xj);
        else
            p = NaN;
        end
        p_values(i,j) = p;
        p_values(j,i) = p;
    end
end
end
