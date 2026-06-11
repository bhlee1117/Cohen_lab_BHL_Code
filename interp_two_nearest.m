function v_new = interp_two_nearest(coords, values, query_coords)
% interp_two_nearest  Linear interpolation from the two closest source points.
%
% For each query coordinate, finds the two nearest original coordinates and
% linearly interpolates: the closer point gets proportionally more weight.
%
%   weight_1 = d2 / (d1 + d2)
%   weight_2 = d1 / (d1 + d2)
%   v_new    = weight_1 * v1 + weight_2 * v2
%
% where d1 <= d2 are the distances to the first and second nearest points.
%
% Usage:
%   v_new = interp_two_nearest(coords, values, query_coords)
%
% Inputs:
%   coords        [N x 2]  (x, y) coordinates of original data points
%   values        [N x 1]  scalar value at each original coordinate
%   query_coords  [M x 2]  (x, y) coordinates to interpolate at
%
% Output:
%   v_new         [M x 1]  interpolated values at query_coords

% -------------------------------------------------------------------------
assert(size(coords, 2) == 2,       'coords must be N x 2');
assert(size(query_coords, 2) == 2, 'query_coords must be M x 2');
assert(numel(values) == size(coords, 1), ...
    'values must have one entry per row of coords');
assert(size(coords, 1) >= 2, 'Need at least 2 source points');

values = values(:);   % ensure column vector

% -------------------------------------------------------------------------
% Pairwise distance matrix  [M x N]
%   Each row i holds the distances from query point i to all N source points.
% -------------------------------------------------------------------------
diff_x = query_coords(:,1) - coords(:,1)';   % [M x N]  broadcast
diff_y = query_coords(:,2) - coords(:,2)';   % [M x N]
D      = sqrt(diff_x.^2 + diff_y.^2);        % [M x N]

% -------------------------------------------------------------------------
% Find the two nearest source points for every query point
% -------------------------------------------------------------------------
[D_sorted, idx_sorted] = sort(D, 2);          % sort each row ascending

d1 = D_sorted(:, 1);     % distance to nearest        [M x 1]
d2 = D_sorted(:, 2);     % distance to 2nd nearest    [M x 1]
v1 = values(idx_sorted(:, 1));               % value at nearest
v2 = values(idx_sorted(:, 2));               % value at 2nd nearest

% -------------------------------------------------------------------------
% Inverse-distance linear weights
%   When the query falls exactly on a source point (d1 == 0),
%   return that value directly to avoid division by zero.
% -------------------------------------------------------------------------
denom  = d1 + d2;
w1     = d2 ./ denom;    % closer point gets higher weight
w2     = d1 ./ denom;

v_new          = w1 .* v1 + w2 .* v2;
v_new(d1 == 0) = v1(d1 == 0);               % exact-match override

end
