function drawPValueLines(p_values, yOffset, varargin)
% drawPValueLinesSeparate Draws separated significance lines avoiding overlap.
%
%   drawPValueLinesSeparate(p_values, yOffset)
%       p_values : NxN symmetric matrix of p-values.
%       yOffset  : scalar y-offset.
%
% Optional name-value pairs:
%       'Format'      : 'star' (default) or 'text'
%       'TextYOffset' : height above line for text
%       'FontSize'    : font size
%       'Alpha'       : significance threshold (default: 0.05)
%       'StepHeight'  : vertical spacing between stacked bars

    % --- Parse optional inputs ---
    p = inputParser;
    addParameter(p, 'Format', 'star', @(x) ischar(x) || isstring(x));
    addParameter(p, 'TextYOffset', 0.1, @isnumeric);
    addParameter(p, 'FontSize', 12, @isnumeric);
    addParameter(p, 'Alpha', 0.05, @isnumeric);
    addParameter(p, 'StepHeight', 0.1, @isnumeric);
    parse(p, varargin{:});
    opts = p.Results;

    numGroups = size(p_values, 1);

    % --- Set base height ---
    ylims = ylim;
    if isscalar(yOffset)
        baseHeight = ylims(2) + yOffset;
    else
        error('yOffset must be a scalar.');
    end

    hold on;

    % Initialize stack levels
    levelMatrix = [];  % [x1, x2, level]
    stepY = opts.StepHeight;

    % Loop over all pairs
    for i = 1:numGroups
        for j = i+1:numGroups
            pVal = p_values(i, j);
            if isnan(pVal) || pVal > opts.Alpha
                continue;
            end

            x1 = i;
            x2 = j;

            % --- Find available level ---
            currentLevel = 0;
            while true
                overlaps = false;
                for k = 1:size(levelMatrix, 1)
                    existing = levelMatrix(k,:);
                    % If horizontal ranges overlap and same level, cannot use
                    if ~(x2 < existing(1) || x1 > existing(2)) && existing(3) == currentLevel
                        overlaps = true;
                        break;
                    end
                end
                if overlaps
                    currentLevel = currentLevel + 1;
                else
                    break;
                end
            end

            % Register this new bar
            levelMatrix = [levelMatrix; x1, x2, currentLevel];

            % Calculate Y position
            y = baseHeight + currentLevel * stepY;
            dy = stepY/3;

            % Draw line
            plot([x1 x1 x2 x2], [y y+dy y+dy y], 'k', 'LineWidth', 1);

            % Generate label
            if strcmpi(opts.Format, 'star')
                if pVal < 0.001
                    label = '***';
                elseif pVal < 0.01
                    label = '**';
                elseif pVal < 0.05
                    label = '*';
                else
                    label = 'n.s.';
                end
            else
                label = sprintf('p=%.3g', pVal);
            end

            % Place label
            text(mean([x1 x2]), y + opts.TextYOffset, label, ...
                'HorizontalAlignment', 'center', ...
                'FontSize', opts.FontSize, ...
                'FontWeight', 'bold');
        end
    end

    hold off;
end