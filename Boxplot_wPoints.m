function p=Boxplot_wPoints(vector1, vector2, faceColor, varargin)
% Boxplot_wPoints Plots box plots with overlaid data points and displays a t-test.
%
%   Boxplot_wPoints(vector1, vector2, faceColor) creates a plot where
%   vector1 and vector2 are plotted as two groups with box charts using the
%   specified faceColor. The function overlays individual data points with a
%   small horizontal jitter and performs an independent t-test (ttest2).
%
%   Boxplot_wPoints(vector1, vector2, faceColor, 'pairwise') runs a paired t-test (ttest)
%   and connects paired data points with lines.
%
%   Inputs:
%       vector1   - Numeric vector for group 1.
%       vector2   - Numeric vector for group 2.
%       faceColor - Color specification (single color or 2×3 RGB matrix).
%       'pairwise' - (Optional) Uses paired t-test (ttest) instead of independent t-test (ttest2).
%
%   Example:
%       v1 = randn(30,1);
%       v2 = v1 + 0.5;  % Paired data
%       Boxplot_wPoints(v1, v2, 'b', 'pairwise');

    % ---------------------------------------------------------------------
    %% 1) Check if 'pairwise' is passed
    % ---------------------------------------------------------------------
    isPairwise = any(strcmp(varargin, 'pairwise'));

    % ---------------------------------------------------------------------
    %% 2) Convert faceColor to RGB if it's a string
    % ---------------------------------------------------------------------
    if ischar(faceColor) || isstring(faceColor)
        faceColor = convertColor(faceColor);
    end
    
    % Ensure faceColor is a 2×3 RGB matrix
    if size(faceColor,1) == 1
        faceColor = [faceColor; faceColor]; % Duplicate if only one color is given
    elseif size(faceColor,1) ~= 2
        error('faceColor must be either a single 1×3 row or a 2×3 matrix.');
    end

    % ---------------------------------------------------------------------
    %% 3) Create Boxplot
    % ---------------------------------------------------------------------        
hold all;
    if exist('boxchart', 'file') % MATLAB R2020b+
        % Create box charts for each group
        boxchart(ones(size(vector1)), vector1, 'BoxFaceColor', faceColor(1,:), 'MarkerStyle', 'none');
        boxchart(2*ones(size(vector2)), vector2, 'BoxFaceColor', faceColor(2,:), 'MarkerStyle', 'none');
    else % Older MATLAB versions
        data = [vector1(:); vector2(:)];
        groups = [repmat(1, length(vector1), 1); repmat(2, length(vector2), 1)];
        boxplot(data, groups, 'Colors', 'k', 'Symbol', '');
    end

    % ---------------------------------------------------------------------
    %% 4) Overlay scatter points with jitter
    % ---------------------------------------------------------------------
    jitterAmount = 0.1;
    
    % Scatter plot for group 1
    x1 = ones(size(vector1)) + jitterAmount * (rand(size(vector1)) - 0.5);
    scatter(x1, vector1, 40, faceColor(1,:), 'filled', 'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', 'k');

    % Scatter plot for group 2
    x2 = 2*ones(size(vector2)) + jitterAmount * (rand(size(vector2)) - 0.5);
    scatter(x2, vector2, 40, faceColor(2,:), 'filled', 'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', 'k');

    % ---------------------------------------------------------------------
    %% 5) Connect paired points with lines (if 'pairwise' is specified)
    % ---------------------------------------------------------------------
    if isPairwise
        if length(vector1) ~= length(vector2)
            error('For a paired t-test, vector1 and vector2 must have the same length.');
        end
        
        for i = 1:length(vector1)
            plot([x1(i), x2(i)], [vector1(i), vector2(i)], '-','color',[0.7 0.7 0.7], 'LineWidth', 1);
        end
    end

    % ---------------------------------------------------------------------
    %% 6) Perform Statistical Test (ttest2 or ttest)
    % ---------------------------------------------------------------------
    if ~isempty(vector2)
        if isPairwise
            [~, p] = ttest(vector1, vector2);
            testType = 'Paired t-test';
        else
            [~, p] = ttest2(vector1, vector2);
            p = ranksum(vector1, vector2);
            %testType = 'Two-sample t-test';
            testType = 'Rank-sum';
        end
        
        % Display p-value
        yLimits = get(gca, 'YLim');
        text(1.5, yLimits(2) - 0.05*range(yLimits), ...
             sprintf('%s: p = %.3f', testType, p), ...
             'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
    end
    
    % ---------------------------------------------------------------------
    %% 7) Format the plot
    % ---------------------------------------------------------------------
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Group 1', 'Group 2'});
    ylabel('Value');
    hold off;
end

%% Helper function: Convert color names to RGB
function rgb = convertColor(colorName)
    switch lower(colorName)
        case {'red','r'}
            rgb = [1 0 0];
        case {'green','g'}
            rgb = [0 1 0];
        case {'blue','b'}
            rgb = [0 0 1];
        case {'cyan','c'}
            rgb = [0 1 1];
        case {'magenta','m'}
            rgb = [1 0 1];
        case {'yellow','y'}
            rgb = [1 1 0];
        case {'black','k'}
            rgb = [0 0 0];
        case {'white','w'}
            rgb = [1 1 1];
        otherwise
            error('Unknown color specification. Use an RGB triplet or a standard color name.');
    end
end
