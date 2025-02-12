function filteredData = pcafilterTrace(data, numComponents)
% PCAFILTEROMITNAN Filters an N x T matrix using PCA, omitting frames with NaN values.
%
% filteredData = pcaFilterOmitNaN(data, numComponents)
%
% INPUT:
%   data          - N x T input matrix (N features, T time points), may contain NaN values.
%   numComponents - Number of principal components to retain.
%
% OUTPUT:
%   filteredData  - PCA-filtered data matrix (N x T), with NaN frames restored.

    % Find frames (columns) with NaN values
    validFrames = all(~isnan(data), 1); % Logical index of valid (non-NaN) frames

    % Omit frames with NaN values
    dataValid = data(:, validFrames);

    % Center the data by subtracting the mean
    dataMean = mean(dataValid, 2); % Mean along time (T)
    centeredData = dataValid - dataMean;

    % Perform singular value decomposition (SVD)
    [U, S, V] = svd(centeredData, 'econ');

    % Retain the specified number of components
    U_reduced = U(:, 1:numComponents); % Leading components
    S_reduced = S(1:numComponents, 1:numComponents); % Reduced singular values
    V_reduced = V(:, 1:numComponents); % Right singular vectors

    % Reconstruct the data using retained components
    filteredDataValid = U_reduced * S_reduced * V_reduced';

    % Add back the mean
    filteredDataValid = filteredDataValid + dataMean;

    % Restore filtered data to the original size, keeping NaN frames
    filteredData = NaN(size(data)); % Initialize output with NaN
    filteredData(:, validFrames) = filteredDataValid; % Fill valid frames
end