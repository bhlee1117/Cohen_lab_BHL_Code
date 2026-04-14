function h = drawScaleBar(barLength, direction, varargin)
% drawScaleBar - Draws a simple scale bar on the current axes.

% h = drawScaleBar(barLength, direction)
% 
% Inputs:
%   barLength : Length of the bar in data units (e.g., pixels)
%   direction : 'horizontal' or 'vertical'
% 
% Optional name-value pairs:
%   'Position' : [x, y] starting coordinate (default: bottom-left inset)
%   'Color'    : Line color (default: 'k')
%   'LineWidth': Width of the line (default: 2)

% Parse inputs
p = inputParser;
addParameter(p, 'Position', []);
addParameter(p, 'Color', 'w');
addParameter(p, 'LineWidth', 2);
parse(p, varargin{:});
opts = p.Results;

ax = gca;
xlim_ = xlim(ax);
ylim_ = ylim(ax);

% Default position: 5% inset from bottom-right
if isempty(opts.Position)
    %x0 = xlim_(1) + 0.75 * range(xlim_);
    x1 = xlim_(1) + 0.95 * range(xlim_);
    y0 = ylim_(1) + 0.9 * range(ylim_);
else
    x1 = opts.Position(1);
    y0 = opts.Position(2);
end

% Define endpoint
switch lower(direction)
    case 'horizontal'
        x0 = x1-barLength;
        y1 = y0;
    case 'vertical'
        x0 = x1;
        y1 = y0 + barLength;
    otherwise
        error('Direction must be ''horizontal'' or ''vertical''.');
end

% Draw the bar
h = line([x0, x1], [y0, y1], ...
    'Color', opts.Color, ...
    'LineWidth', opts.LineWidth);
end
