function coloredImage = showScaleImage(binaryStack, valueVector, colormapName, valueRange)
% colorScaleImageStack - Applies color scaling to a binary image stack based on a vector of values and outputs a single colored image.
%
% Syntax: coloredImage = colorScaleImageStack(binaryStack, valueVector, colormapName, valueRange)
%
% Inputs:
%   binaryStack - 3D binary image stack (X x Y x N)
%   valueVector - 1xN vector with corresponding values for each image slice
%   colormapName - String specifying the colormap (e.g., 'jet', 'parula', etc.)
%   valueRange - 1x2 vector specifying the [min, max] range of values for the color axis (optional)
%
% Outputs:
%   coloredImage - 3D array (X x Y x 3) representing the final combined colorized image

% Validate inputs
if size(binaryStack, 3) ~= length(valueVector)
    error('The number of slices in binaryStack must match the length of valueVector.');
end

% Get the specified colormap
if nargin < 3 || isempty(colormapName)
    colormapName = 'turbo'; % Default colormap
end
cmap = colormap(colormapName);

% Determine value range for normalization
if nargin < 4 || isempty(valueRange)
    valueMin = min(valueVector);
    valueMax = max(valueVector);
else
    valueMin = valueRange(1);
    valueMax = valueRange(2);
end

% Normalize the valueVector to map to colormap indices
normalizedValues = (valueVector - valueMin) / (valueMax - valueMin);
normalizedValues = max(0, min(1, normalizedValues)); % Clamp values to [0, 1]
colormapIndices = round(normalizedValues * (size(cmap, 1) - 1)) + 1;

% Preallocate the combined colored image
[X, Y, ~] = size(binaryStack);
coloredImage = zeros(X, Y, 3);

% Keep track of overlapping regions
combinedMask = zeros(X, Y);

% Apply colormap and combine slices into a single image
for i = 1:size(binaryStack, 3)
    % Extract the binary slice and its corresponding color
    binarySlice = binaryStack(:, :, i);
    binarySlice(combinedMask > 0) = 0; % Omit overlapping regions
    color = cmap(colormapIndices(i), :); % RGB color from colormap
    
    % Add the weighted color slice to the final image
    for channel = 1:3
        coloredImage(:, :, channel) = coloredImage(:, :, channel) + binarySlice * color(channel);
    end
    
    % Update the combined mask
    combinedMask = combinedMask + binarySlice;
end

% Normalize the combined image to avoid over-saturation
coloredImage = coloredImage / max(coloredImage(:));

% Display the final colored image to verify output
imagesc(coloredImage);
end
