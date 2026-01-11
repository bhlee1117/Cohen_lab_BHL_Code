function videoOut = MovMaskNaN(videoIn, mask)
% maskVideoWithNaN - Sets pixel values to NaN where mask == 1 across all frames
%
%   videoOut = MovMaskNaN(videoIn, mask)
%
%   Inputs:
%       videoIn : [X x Y x T] numeric array (e.g., fluorescence video)
%       mask    : [X x Y] logical or binary matrix (1 = pixels to set to NaN)
%
%   Output:
%       videoOut : [X x Y x T] same as input, but masked pixels are set to NaN

    % Check size compatibility
    if ~isequal(size(videoIn,1), size(mask,1)) || ~isequal(size(videoIn,2), size(mask,2))
        error('The spatial dimensions of videoIn and mask must match.');
    end

    % Convert mask to logical if not already
    mask = logical(mask);

    % Replicate mask across time dimension
    mask3D = repmat(mask, [1 1 size(videoIn,3)]);

    % Copy video and set masked pixels to NaN
    videoOut = videoIn;
    videoOut(mask3D) = NaN;

end
