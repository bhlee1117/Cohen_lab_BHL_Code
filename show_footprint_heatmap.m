function colored_image = show_footprint_heatmap(c_ftprnt, dat, dat_min, dat_max)
% colored_image = show_footprint_heatmap(c_ftprnt, dat, dat_min, dat_max)
% show_footprint_heatmap  Overlay cell footprints colored by a scalar value.
%
% Renders a false-color RGB image in which every cell (ROI) is filled with
% the color that corresponds to its scalar value, mapped through the
% 'turbo' colormap.  The result is displayed in the current tile / axes
% and also returned as an [nY x nX x 3] RGB image.
%
% Syntax:
%   colored_image = show_footprint_heatmap(c_ftprnt, dat)
%   colored_image = show_footprint_heatmap(c_ftprnt, dat, dat_min, dat_max)
%
% Inputs:
%   c_ftprnt  [nY x nX x nCells]  cell footprint stack (logical or numeric).
%             Any nonzero pixel is treated as belonging to the cell.
%   dat       [nCells x 1] or [1 x nCells]  scalar value per cell
%             (e.g. Fslope, dF/F, correlation coefficient).
%   dat_min   (optional) scalar — colormap lower bound. Default: min(dat).
%   dat_max   (optional) scalar — colormap upper bound. Default: max(dat).
%
% Output:
%   colored_image  [nY x nX x 3]  RGB image. Background pixels are black (0).
%                  Overlapping cell footprints have their colors summed.
%
% Dependencies (Cohen Lab utilities):
%   grs2rgb  — maps a value vector to RGB colors given a colormap and range
%   imshow2  — displays an image in the current axes
%
% Example:
%   figure;
%   tiledlayout(2, 1);
%   nexttile;
%   img = show_footprint_heatmap(ftprnt, slopes);
%   nexttile;
%   img = show_footprint_heatmap(ftprnt, Rsq_values, 0, 1);

% -------------------------------------------------------------------------
% Input validation
% -------------------------------------------------------------------------
if isempty(c_ftprnt) || isempty(dat)
    error('show_footprint_heatmap:emptyInput', ...
        'c_ftprnt and dat must not be empty.');
end

if ndims(c_ftprnt) ~= 3
    error('show_footprint_heatmap:badInput', ...
        'c_ftprnt must be a 3-D array [nY x nX x nCells], but got %d-D.', ...
        ndims(c_ftprnt));
end

if ~isvector(dat)
    error('show_footprint_heatmap:badInput', ...
        'dat must be a vector [nCells x 1], but got size [%s].', ...
        num2str(size(dat)));
end

nCells = size(c_ftprnt, 3);
if nCells ~= numel(dat)
    error('show_footprint_heatmap:sizeMismatch', ...
        'c_ftprnt has %d cells (3rd dim) but dat has %d elements.', ...
        nCells, numel(dat));
end

dat = dat(:);   % enforce column vector

% --- color range bounds ---------------------------------------------------
if nargin < 3 || isempty(dat_min)
    dat_min = min(dat);
end
if nargin < 4 || isempty(dat_max)
    dat_max = max(dat);
end

if ~isscalar(dat_min) || ~isscalar(dat_max)
    error('show_footprint_heatmap:badInput', ...
        'dat_min and dat_max must be scalars.');
end

if dat_min > dat_max
    error('show_footprint_heatmap:badRange', ...
        'dat_min (%.4g) must be <= dat_max (%.4g).', dat_min, dat_max);
end

if dat_min == dat_max
    warning('show_footprint_heatmap:flatRange', ...
        'dat_min equals dat_max (%.4g): all cells will share the same color.', ...
        dat_min);
    dat_max = dat_min + 1;   % avoid division by zero inside grs2rgb
end

% -------------------------------------------------------------------------
% Build colored image
% -------------------------------------------------------------------------
% Binarize footprints: any nonzero pixel counts as part of the cell
c_ftprnt = double(c_ftprnt > 0);   % [nY x nX x nCells]

% Map each cell's scalar value to an RGB color via the turbo colormap.
% grs2rgb returns a [nCells x 3] matrix after squeeze.
colr = squeeze(grs2rgb(dat, colormap(turbo), dat_min, dat_max));  % [nCells x 3]

% Combine: broadcast footprint [nY x nX x nCells x 1] against
% color     [1  x  1  x nCells x 3], sum over the cell dimension.
% Result after squeeze: [nY x nX x 3] RGB image.
colored_image = squeeze(sum(c_ftprnt .* reshape(colr, 1, 1, nCells, 3), 3));

% -------------------------------------------------------------------------
% Display
% -------------------------------------------------------------------------
ax = nexttile;
imshow2(colored_image, []);
title(ax, sprintf('Footprint heatmap  (%d cells,  range [%.3g, %.3g])', ...
    nCells, dat_min, dat_max), 'FontSize', 9);

end
