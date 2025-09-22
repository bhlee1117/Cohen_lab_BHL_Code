function drawPValueLines(p_values, yOffset, varargin)
% drawPValueLines Draws significance lines for p-value matrix.
%   drawPValueLines(p_values, yOffset, 'XCoord', x)
%
%   - p_values: NxN symmetric matrix
%   - yOffset: scalar for vertical offset above current axis limit
%
% Optional name-value pairs:
%   'XCoord'      : 1xN vector of x-axis coordinates
%   'Format'      : 'star' or 'text'
%   'TextYOffset' : vertical offset for label
%   'FontSize'    : font size for label
%   'Alpha'       : significance threshold
%   'StepHeight'  : vertical spacing between stacked lines

    % --- Parse optional inputs ---
    p = inputParser;
    addParameter(p, 'XCoord', [], @(x) isnumeric(x) && isvector(x));
    addParameter(p, 'Format', 'star', @(x) ischar(x) || isstring(x));
    addParameter(p, 'TextYOffset', 0.1, @isnumeric);
    addParameter(p, 'FontSize', 12, @isnumeric);
    addParameter(p, 'Alpha', 0.05, @isnumeric);
    addParameter(p, 'StepHeight', 0.1, @isnumeric);
    parse(p, varargin{:});
    opts = p.Results;

    numGroups = size(p_values, 1);

    % --- Set default X coordinates if not provided ---
    if isempty(opts.XCoord)
        xCoord = 1:numGroups;
    else
        xCoord = opts.XCoord;
        if numel(xCoord) ~= numGroups
            error('Length of XCoord must match size of p-value matrix.');
        end
    end

    % --- Set base height ---
    ylims = ylim;
    if isscalar(yOffset)
        baseHeight = ylims(2) + yOffset;
    else
        error('yOffset must be a scalar.');
    end

    hold on;
    levelMatrix = [];  % [x1, x2, level]
    stepY = opts.StepHeight;

    % Loop over all pairs
    for i = 1:numGroups
        for j = i+1:numGroups
            pVal = p_values(i, j);
            if isnan(pVal) || pVal > opts.Alpha
                continue;
            end

            x1 = xCoord(i);
            x2 = xCoord(j);

            % --- Find available level ---
            currentLevel = 0;
            while true
                overlaps = false;
                for k = 1:size(levelMatrix, 1)
                    existing = levelMatrix(k,:);
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

            % Register this line
            levelMatrix = [levelMatrix; x1, x2, currentLevel];

            % --- Draw line ---
            y = baseHeight + currentLevel * stepY;
            dy = stepY/3;
            plot([x1 x1 x2 x2], [y y+dy y+dy y], 'k', 'LineWidth', 1);

            % --- Label ---
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

            text(mean([x1 x2]), y + opts.TextYOffset, label, ...
                'HorizontalAlignment', 'center', ...
                'FontSize', opts.FontSize, ...
                'FontWeight', 'bold');
        end
    end

    hold off;
end
