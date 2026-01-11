function [density, x_edges, y_edges] = scatter_heatmap(x, y, Nbinx,Nbiny)
% SCATTER_HEATMAP Creates a 2D heatmap from scattered x-y data.
%
%   [counts, x_edges, y_edges] = SCATTER_HEATMAP(x, y, Nbin) computes and
%   displays a heatmap representing the density of (x, y) data points.
%
%   INPUTS:
%     x      - Vector of x-coordinates (must be same length as y)
%     y      - Vector of y-coordinates
%     Nbin   - Number of bins to use along each axis (scalar)
%
%   OUTPUTS:
%     counts    - 2D matrix of normalized point counts per bin
%     x_edges   - Bin edges along the x-axis
%     y_edges   - Bin edges along the y-axis
%
%   DESCRIPTION:
%     This function partitions the x-y space into a 2D grid of Nbin × Nbin
%     equally spaced bins, computes the number of points falling into each
%     bin (normalized by total number of points), and displays the result
%     as a heatmap. The density is color-coded and visualized using `imagesc`.
%
%   NOTES:
%     - The heatmap is scaled such that the sum of all bins equals 1.
%     - The colorbar indicates the relative density in each bin.
%
%   EXAMPLE:
%     x = randn(1000,1);
%     y = 0.5*x + randn(1000,1);
%     scatter_heatmap(x, y, 50);
%
%   See also: HISTCOUNTS2, IMAGESC, COLORBAR

% Define the bin edges
% x_edges = linspace(min(x), max(x), Nbin);
% y_edges = linspace(min(y), max(y), Nbin);

% Compute the 2D histogram (heatmap data)
% Npoints = length(x);
% counts = histcounts2(x, y, x_edges, y_edges) / Npoints;

data = [x(:), y(:)];
xc=linspace(min(x), max(x), Nbinx);
yc=linspace(min(y), max(y), Nbiny);
[xx, yy] = ndgrid(linspace(min(x), max(x), Nbinx), ...
    linspace(min(y), max(y), Nbiny));
[f,~] = ksdensity(data, [xx(:), yy(:)]);
%Fgrid = griddata(xx(:), yy(:), f, data(:,1), data(:,2), 'nearest');
density = reshape(f,Nbinx,Nbiny);

% Plot the heatmap
% imagesc(x_edges, y_edges, counts');
imagesc(xc, yc, density');
axis xy;
colorbar;
end
