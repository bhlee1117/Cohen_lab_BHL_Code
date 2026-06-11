function [coord1D principal_axis]= projectTrunkaxis(X)
% PROJECT_TO_NEURONAL_AXIS - Projects 2D subcellular coordinates onto a 
% 1D principal axis, with soma at 0, apical dendrite positive, basal negative.
%
% INPUT:
%   X : N x 2 matrix of [x, y] coordinates, with X(1,:) as the soma.
%
% OUTPUT:
%   coord1D : N x 1 vector of 1D coordinates along the apical-basal axis

    if size(X, 2) ~= 2
        error('Input must be an N x 2 matrix.');
    end

    soma = X(1, :);
    X_centered = X - soma;  % Center data at soma

    % PCA to get the principal axis
    [coeff, ~, ~] = pca(X_centered);
    principal_axis = coeff(:,1);  % First principal component

    % Project onto principal axis
    projections = X_centered * principal_axis;

    % Determine correct polarity (apical is further than basal)
    maxProj = max(projections);
    minProj = min(projections);
    if abs(minProj) > abs(maxProj)
        projections = -projections;
    end

    % Soma should be zero
    coord1D = projections;
end
