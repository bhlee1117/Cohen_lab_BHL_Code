function [xx, yy] = get_figsize(fig)
% GET_FIGSIZE  Get the current (or specified) figure size in millimeters.
%
% Returns the figure's on-screen size converted to millimeters, using the
% same screen-DPI conversion that SET_FIGSIZE uses to go from millimeters
% to pixels. This makes the two functions round-trip compatible:
%
%   [xx, yy] = get_figsize(fig);
%   set_figsize(xx, yy, fig);   % restores the same on-screen size
%
% Usage:
%   [xx, yy] = get_figsize()      - Read from the current figure
%   [xx, yy] = get_figsize(fig)   - Read from a specified figure handle
%   sz       = get_figsize(...)   - Single output: sz = [xx yy]
%
% Inputs:
%   fig - (optional) Figure handle. Defaults to gcf.
%
% Outputs:
%   xx  - Figure width  in millimeters
%   yy  - Figure height in millimeters

if nargin < 1 || isempty(fig)
    fig = gcf;
end

% Figure position in pixels: [left bottom width height]
pos = get(fig, 'Position');

% Same screen DPI conversion used by set_figsize
dpi   = get(0, 'ScreenPixelsPerInch');
mm2px = dpi / 25.4;

xx = pos(3) / mm2px;
yy = pos(4) / mm2px;

% Allow sz = get_figsize(...) single-output usage
if nargout <= 1
    xx = [xx, yy];
end

end
