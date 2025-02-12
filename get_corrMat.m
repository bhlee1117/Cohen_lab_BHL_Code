function rowCorrelations = get_corrMat(matrix1, matrix2, frameRange)
    % calculateRowCorrelations calculates correlation values between rows of two matrices.
    %
    % Inputs:
    %   - matrix1: An NxT matrix with potential NaN values.
    %   - matrix2: A KxT matrix with potential NaN values.
    %   - frameRange: A vector specifying the columns (frames) to include in the calculation.
    %
    % Output:
    %   - rowCorrelations: An NxK matrix containing pairwise row correlations.

    % Check inputs
    if nargin < 3 || isempty(frameRange)
        frameRange = 1:size(matrix1, 2); % Use all frames if none are specified
    end

    % Extract the specified frames
    matrix1Subset = matrix1(:, frameRange);
    matrix2Subset = matrix2(:, frameRange);

    % Initialize the correlation matrix
    N = size(matrix1, 1);
    K = size(matrix2, 1);
    rowCorrelations = nan(N, K);

    % Compute pairwise correlations
    for i = 1:N
        for j = 1:K
            % Get rows i and j
            row1 = matrix1Subset(i, :);
            row2 = matrix2Subset(j, :);

            % Remove NaN values
            validIndices = ~isnan(row1) & ~isnan(row2);
            validRow1 = row1(validIndices);
            validRow2 = row2(validIndices);

            % Calculate correlation if enough data points exist
            if length(validRow1) > 1 % At least 2 points needed for correlation
                rowCorrelations(i, j) = corr(validRow1(:), validRow2(:));
            end
        end
    end
end

