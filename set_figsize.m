function set_figsize(xx, yy, fig)
% SET_FIGSIZE  Set the current (or specified) figure size in millimeters and
%              whiten the background for clean export to Adobe Illustrator.
%
% Usage:
%   set_figsize(xx, yy)        - Apply to the current figure
%   set_figsize(xx, yy, fig)   - Apply to a specified figure handle
%
% Inputs:
%   xx  - Figure width  in millimeters
%   yy  - Figure height in millimeters
%   fig - (optional) Figure handle. Defaults to gcf.

if nargin < 3
    fig = gcf;
end

% Convert mm to points (1 mm = 2.8346 points), used by MATLAB's PaperSize
mm2pt = 2.8346;

% Get current figure position to preserve the on-screen top-left corner
pos = get(fig, 'Position');  % [left bottom width height] in pixels

% Screen DPI for pixel conversion (MATLAB default is 96 dpi on most systems)
dpi = get(0, 'ScreenPixelsPerInch');
mm2px = dpi / 25.4;

% Set figure size in pixels (preserves position on screen)
set(fig, 'Position', [pos(1), pos(2), xx * mm2px, yy * mm2px]);

% Set paper size and position for export (in centimeters)
set(fig, 'PaperUnits', 'centimeters');
set(fig, 'PaperSize',     [xx/10, yy/10]);
set(fig, 'PaperPosition', [0, 0, xx/10, yy/10]);
set(fig, 'PaperPositionMode', 'manual');

% Apply light theme for clean Illustrator export (compatible with all versions)
set(fig, 'Color', [1 1 1]);                          % white figure background
set(gca, 'Color', [1 1 1]);                          % white axes background
set(gca, 'XColor', [0 0 0], 'YColor', [0 0 0]);     % black axis lines & labels
set(findall(fig, '-property', 'Color', 'Type', 'text'), 'Color', [0 0 0]);  % all text black
set(gca, 'GridColor', [0.15 0.15 0.15]);             % dark grid lines
set(fig, 'InvertHardcopy', 'off');  % prevent MATLAB from resetting bg on print
copygraphics(gcf,'Resolution',300,'ContentType','image');
end
