function colors = vec2cmap(vec, cmapName, range)
% vec2cmap Maps a numeric vector to RGB colors using a colormap.
%
%   colors = map_to_colormap(vec, cmapName)
%   colors = map_to_colormap(vec, cmapName, range)
%
%   Inputs:
%       vec      - Numeric vector to map (1D array)
%       cmapName - Name of the colormap (e.g., 'turbo', 'parula', 'jet')
%       range    - [min max] range for mapping (optional)
%
%   Output:
%       colors   - N x 3 RGB matrix, one row per value in vec
if nargin<3
    range=[min(vec) max(vec)];
end

    % Get colormap
    N = 256;
    if ischar(cmapName)
    cmap = feval(cmapName, N);
    else
    cmap = cmapName;    
    end

    % Clamp range
    vecClamped = max(min(vec, range(2)), range(1));

    % Normalize to [0 1]
    normVals = (vecClamped - range(1)) / (range(2) - range(1));

    % Map to color index
    idx = max(1, round(normVals * (N-1)) + 1);

    % Map to RGB
    colors = cmap(idx, :);
end
