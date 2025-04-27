function contour_cmap(X, Y, Z, n_levels, color_list)
% contour_with_colors plots contour lines with specified line colors
%
% Inputs:
%   X, Y, Z       - meshgrid and scalar field
%   n_levels      - number of contour levels
%   color_list    - optional:
%                   (a) N x 3 RGB matrix
%                   (b) cell array of color chars like {'r','g','b'}
%                   If omitted, uses 'jet' colormap.

% Compute contour levels
zmin = min(Z(:));
zmax = max(Z(:));
levels = linspace(zmin, zmax, n_levels);

% Determine colors
if nargin < 5
    color_list = jet(n_levels);  % default colormap
end

% Plot each level separately
hold on;
for i = 1:n_levels
    if iscell(color_list)
        c = color_list{i};
    else
        c = color_list(i, :);
    end
    contour(X, Y, Z, [levels(i) levels(i)], 'LineColor', c, 'LineWidth', 1.2);
end
hold off;
end