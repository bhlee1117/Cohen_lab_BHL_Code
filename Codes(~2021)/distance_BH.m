function D = distance_BH(X, Y)
    % COMPUTE_DISTANCE computes pairwise Euclidean distances between points in X and Y.
    %
    %   D = compute_distance(X, Y)
    %
    %   Inputs:
    %       X - N x D matrix, where each row is a coordinate (2D: (x, y) or 3D: (x, y, z)).
    %       Y - M x D matrix, where each row is a coordinate (2D: (x, y) or 3D: (x, y, z)).
    %
    %   Output:
    %       D - N x M matrix, where D(i,j) is the Euclidean distance between X(i,:) and Y(j,:).
    
    % Ensure both matrices have the same number of columns (dimensions)
    if size(X, 2) ~= size(Y, 2)
        error('X and Y must have the same number of columns (2 for (x, y) or 3 for (x, y, z)).');
    end
    
    % Compute Euclidean distance using vectorized operations
    X_exp = permute(X, [1, 3, 2]); % Convert X to N x 1 x D
    Y_exp = permute(Y, [3, 1, 2]); % Convert Y to 1 x M x D
    
    % Compute pairwise distance
    D = sqrt(sum((X_exp - Y_exp) .^ 2, 3));
end
