function filteredData = pcafilterTrace_template(data, templates)
% PCAFILTEROMITNAN Filters an N x T matrix using PCA, omitting frames with NaN values.
%
% filteredData = pcaFilterOmitNaN(data, numComponents)
%
% INPUT:
%   data          - N x T input matrix (N features, T time points), may contain NaN values.
%   components2remain - PC components to take
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

% [ysize, xsize, nframes] = size(mov_res);
% 
% avgimg = mean(mov_res, 3);
% mov2d = reshape(mov_res, [ysize*xsize, nframes]);
% avgimg2d = reshape(avgimg, [ysize*xsize, 1]);
% dmov2d = mov2d - repmat(avgimg2d, [1 nframes]);
% covmat = dmov2d'*dmov2d;
% [V, D] = eig(covmat);
% eigvals = flipud(diag(D));
% V = V(:,end:-1:1);
% vSign = sign(max(V) - max(-V));  % make the largest value always positive
% V = V.*vSign;
% eigvecs = V;
EigV=(templates);

[Q, ~] = qr(EigV, 0); % Orthonormal basis in columns of Q
coeffs = Q'*data;

filteredData = (Q*coeffs);

for n=1:size(coeffs,1)
[fraction(n)] = get_variance(dataValid, coeffs(n,:));
end
disp(['Variance explained by template is ' num2str(sum(fraction))]);