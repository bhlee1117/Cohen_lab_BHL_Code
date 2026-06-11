function [bin_sums, x_edges, y_edges] = binning_data2D(x_coords, y_coords, values, x_edges, y_edges)
% binning_data2D  Accumulate values into 2D spatial bins defined by custom edges.
%
% Assigns each (x, y) data point to a 2D bin and sums the corresponding
% values within each bin. Points outside the bin edge range are ignored.
%
% Usage:
%   [bin_sums, x_edges, y_edges] = binning_data2D(x_coords, y_coords, values, x_edges, y_edges)
%
% Inputs:
%   x_coords  - [N x 1] x-coordinate of each data point
%   y_coords  - [N x 1] y-coordinate of each data point
%   values    - [N x 1] value associated with each data point (e.g. amplitude)
%   x_edges   - [1 x (nXbins+1)] bin edges along the x-axis
%   y_edges   - [1 x (nYbins+1)] bin edges along the y-axis
%
% Outputs:
%   bin_sums  - [nYbins x nXbins] matrix of summed values per 2D bin
%   x_edges   - x bin edges (same as input, returned for convenience)
%   y_edges   - y bin edges (same as input, returned for convenience)
%
% Note:
%   Points outside the range of x_edges or y_edges are silently discarded.
%   Bins are left-inclusive: bin (i,j) spans [x_edges(j), x_edges(j+1)) x
%   [y_edges(i), y_edges(i+1)).
%
% Example:
%   x_edges = linspace(0, 1, 11);   % 10 bins along x
%   y_edges = linspace(0, 1, 11);   % 10 bins along y
%   bin_sums = binning_data2D(x, y, amplitude, x_edges, y_edges);
%   imagesc(bin_sums);
%
% See also: binning_data, binning_data_median

if numel(x_coords) ~= numel(y_coords) || numel(x_coords) ~= numel(values)
    error('binning_data2D: x_coords, y_coords, and values must all have the same number of elements.');
end

nXbins = numel(x_edges) - 1;
nYbins = numel(y_edges) - 1;

% Assign each point to a bin index (NaN = out of range)
x_bin_idx = discretize(x_coords, x_edges);
y_bin_idx = discretize(y_coords, y_edges);

% Keep only points that fall within the bin range on both axes
in_range = ~isnan(x_bin_idx) & ~isnan(y_bin_idx);
x_bin_idx = x_bin_idx(in_range);
y_bin_idx = y_bin_idx(in_range);
values    = values(in_range);

% Accumulate values into 2D bin grid [nYbins x nXbins]
bin_sums = zeros(nYbins, nXbins);
for iPoint = 1:numel(values)
    bin_sums(y_bin_idx(iPoint), x_bin_idx(iPoint)) = ...
        bin_sums(y_bin_idx(iPoint), x_bin_idx(iPoint)) + values(iPoint);
end
end
