function Boxplot_wPoints(vector1, vector2, faceColor)
% Boxplot_wPoints Plots box plots with overlaid data points and displays a t-test.
%
%   Boxplot_wPoints(vector1, vector2, faceColor) creates a plot where
%   vector1 and vector2 are plotted as two groups with box charts using the
%   specified faceColor. The function overlays individual data points with a
%   small horizontal jitter and performs an independent t-test, displaying the
%   p-value on the plot.
%
%   Inputs:
%       vector1   - Numeric vector for group 1.
%       vector2   - Numeric vector for group 2.
%       faceColor - Color specification. Options:
%                     (a) A single color (char/string or 1×3 RGB) that
%                         applies to both groups.
%                     (b) A 2×3 RGB matrix to specify separate colors for
%                         each group, e.g. [R1 G1 B1; R2 G2 B2].
%
%   Example:
%       v1 = randn(30,1);
%       v2 = randn(35,1) + 0.5;
%       % Single color
%       Boxplot_wPoints(v1, v2, 'r');
%       % Two colors
%       Boxplot_wPoints(v1, v2, [0.2 0.6 0.8; 0.7 0.3 0.3]);

    % ---------------------------------------------------------------------
    %% 1) Convert faceColor to an Nx3 RGB matrix if it is a char/string
    % ---------------------------------------------------------------------
    if ischar(faceColor) || isstring(faceColor)
        switch lower(faceColor)
            case {'red','r'}
                faceColor = [1 0 0];
            case {'green','g'}
                faceColor = [0 1 0];
            case {'blue','b'}
                faceColor = [0 0 1];
            case {'cyan','c'}
                faceColor = [0 1 1];
            case {'magenta','m'}
                faceColor = [1 0 1];
            case {'yellow','y'}
                faceColor = [1 1 0];
            case {'black','k'}
                faceColor = [0 0 0];
            case {'white','w'}
                faceColor = [1 1 1];
            otherwise
                error('Unknown color specification. Please use an RGB triplet or a standard color name.');
        end
    end
    
    % ---------------------------------------------------------------------
    %% 2) Handle different shapes of faceColor (1×3 or 2×3)
    % ---------------------------------------------------------------------
    if size(faceColor,1) == 1
        % If only one color is given, repeat for both groups.
        faceColor = [faceColor; faceColor];
    elseif size(faceColor,1) ~= 2
        error(['faceColor must be either a single 1×3 row or ' ...
               'a 2×3 matrix for two groups.']);
    end
    
    % Now faceColor(1,:) is the color for group 1, and
    %     faceColor(2,:) is the color for group 2.

    % ---------------------------------------------------------------------
    %% 3) Check if boxchart is available (R2020b or later)
    % ---------------------------------------------------------------------
    if exist('boxchart', 'file')
        %figure;
        hold on;
        
        % Create box charts for each group
        % For group 1 (position 1)
        bc1 = boxchart(ones(size(vector1)), vector1, ...
                       'BoxFaceColor', faceColor(1,:), ...
                       'MarkerStyle', 'none');
        % For group 2 (position 2)
        bc2 = boxchart(2*ones(size(vector2)), vector2, ...
                       'BoxFaceColor', faceColor(2,:), ...
                       'MarkerStyle', 'none');
                   
    else
        % Fallback to using legacy boxplot
        data = [vector1(:); vector2(:)];
        groups = [repmat(1, length(vector1), 1); repmat(2, length(vector2), 1)];
        
        %figure;
        boxplot(data, groups, 'Colors', 'k', 'Symbol', '');
        hold on;
        % In older MATLAB versions, setting the box fill color is not as straightforward.
        % This code won't color the old 'boxplot' with faceColor directly, but
        % the scatter below will still show the correct face colors.
    end

    % ---------------------------------------------------------------------
    %% 4) Overlay scatter points with horizontal jitter
    % ---------------------------------------------------------------------
    jitterAmount = 0.1; % Adjust if needed
    
    % Group 1 scatter
    x1 = ones(size(vector1)) + jitterAmount*(rand(size(vector1)) - 0.5);
    scatter(x1, vector1, ...
        40, ...                      % marker size
        faceColor(1,:), 'filled', ...% marker face color
        'MarkerFaceAlpha', 0.7, ...
        'MarkerEdgeColor', 'k');

    % Group 2 scatter
    x2 = 2*ones(size(vector2)) + jitterAmount*(rand(size(vector2)) - 0.5);
    scatter(x2, vector2, ...
        40, ...
        faceColor(2,:), 'filled', ...
        'MarkerFaceAlpha', 0.7, ...
        'MarkerEdgeColor', 'k');
    
    % ---------------------------------------------------------------------
    %% 5) Perform statistical test (independent two-sample t-test) & annotate
    % ---------------------------------------------------------------------
    [~, p] = ttest2(vector1, vector2);
    
    yLimits = get(gca, 'YLim');
    text(1.5, yLimits(2) - 0.05*range(yLimits), ...
         sprintf('p = %.3f', p), ...
         'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
    
    % ---------------------------------------------------------------------
    %% 6) Format the plot
    % ---------------------------------------------------------------------
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Group 1', 'Group 2'});
    ylabel('Value');
    hold off;
end
