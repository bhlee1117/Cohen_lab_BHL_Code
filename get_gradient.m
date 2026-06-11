function [grad_mag, grad_dir, grad_vec] = get_gradient(val, coord, r_neighbor)
% function [grad_mag, grad_dir, grad_vec] = get_gradient(val, coord, r_neighbor)
% Inputs
% val:   N x 1  (your value vector)
% coord: N x 2  (x, y coordinates)

% r_neighbor: radius to define local neighbourhood (tune to your data scale)

N = size(coord, 1);
grad_vec = nan(N, 2);
D2 = pdist2(coord, coord, 'squaredeuclidean');
sigma2 = (r_neighbor/2)^2;

for i = 1:N
    nbr = find(D2(i,:) <= r_neighbor^2);
    if numel(nbr) < 3, continue; end

    % Gaussian distance weight ONLY — symmetric high values cancel correctly
    w = exp(-D2(i, nbr)' / (2 * sigma2));
    W = diag(w);

    xy = coord(nbr,:) - coord(i,:);
    A  = [ones(numel(nbr),1), xy];
    coef = (A'*W*A) \ (A'*W*val(nbr));
    grad_vec(i,:) = coef(2:3)';
end

grad_mag = sqrt(sum(grad_vec.^2, 2));
grad_dir = atan2d(grad_vec(:,2), grad_vec(:,1));

% Scale magnitude by local value — suppress arrows in low-signal regions
% without biasing the gradient direction/magnitude estimate
val_norm = (val - min(val)) / (range(val) + eps);
grad_mag_display = grad_mag .* val_norm;
end