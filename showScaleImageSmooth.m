function RGB = showScaleImageSmooth(mask, V, cmapNameOrArray, valueRange, sigma)
% SHOWSCALEIMAGESMOOTH Smoothly interpolates ROI values into RGB image.
%
% Inputs:
%   mask       - [X x Y x N] binary mask for N ROIs
%   V          - [1 x N] values for each ROI
%   cmap       - colormap name (e.g. 'turbo') or [M x 3] colormap
%   valueRange - [min, max] value range for color scaling (optional)
%   sigma      - Gaussian smoothing width (default = 2)
%
% Output:
%   RGB        - [X x Y x 3] smooth color-coded image

    if nargin < 5 || isempty(sigma)
        sigma = 2;
    end

    if nargin < 4 || isempty(valueRange)
        valueMin = min(V(:));
        valueMax = max(V(:));
    else
        valueMin = valueRange(1);
        valueMax = valueRange(2);
    end

    % Initialize
    [X, Y, N] = size(mask);
    V = double(V(:));
    value_map = zeros(X, Y);
    weight_map = zeros(X, Y);

    % Weighted accumulation
    for i = 1:N
        roi = double(mask(:,:,i));
        value_map = value_map + roi * V(i);
        weight_map = weight_map + roi;
    end

    % Smooth both maps
    value_blur = imgaussfilt(value_map, sigma);
    weight_blur = imgaussfilt(weight_map, sigma);

    % Interpolated value field
    interp_val = value_blur ./ max(weight_blur, eps);

    % Normalize using provided valueRange
    interp_norm = (interp_val - valueMin) / (valueMax - valueMin);
    interp_norm = min(max(interp_norm, 0), 1);  % clamp to [0, 1]

    % Load or use given colormap
    if ischar(cmapNameOrArray) || isstring(cmapNameOrArray)
        cmap = feval(char(cmapNameOrArray), 256);
    else
        cmap = cmapNameOrArray;
    end

    % Map to RGB
    idx = round(1 + interp_norm * (size(cmap,1)-1));
    idx = min(max(idx, 1), size(cmap,1));
    RGB = ind2rgb(idx, cmap);

    % Optional: mask boundary
    if nargout == 0 || ~isempty(getenv('SHOW_ROI_MASK'))  % or custom flag
        boundary = max(mask, [], 3) > 0;
        imshow2(RGB .* boundary,[]);
    end
end
