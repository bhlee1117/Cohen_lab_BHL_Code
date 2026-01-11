function masked = maskEdge(matrix, width, C)
% masked = maskEdge(matrix, width, C)
% Input: matrix = 2D or 3D array (rows x cols x T) to mask edges,
%        width = edge boundary width to zero out (default: 5)
% Output: masked = matrix with edge boundaries set to 0

    if nargin < 2
        width = 5;
    end

    if length(size(matrix))<=3
    masked = matrix;
    masked(1:width, :,:) = C;
    masked(end-width+1:end, :,:) = C;
    masked(:, 1:width,:) = C;
    masked(:, end-width+1:end,:) = C;
    end

    if length(size(matrix))==4
    masked = matrix;
    masked(1:width, :,:,:) = C;
    masked(end-width+1:end, :,:,:) = C;
    masked(:, 1:width,:,:) = C;
    masked(:, end-width+1:end,:,:) = C;
    end
end