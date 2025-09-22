function maskSequence = line2masksequence(lineImage, masks, startPt)
% line2masksequence: Returns sequence of masks that a line contour passes through
%
% Inputs:
%   - lineImage: s1 x s2 binary image with a line (nonzero pixels trace the path)
%   - masks:     s1 x s2 x N logical mask stack
%   - startPt:   [x, y] approximate starting point for line tracing
%   - endPt:     [x, y] (optional) approximate ending point (not used directly)
%
% Output:
%   - maskSequence: vector of mask indices in order encountered along the traced line

    if ~isequal(size(lineImage), size(masks(:,:,1)))
        error('lineImage and masks must have same spatial dimensions.');
    end

    % Trace the line from the binary image using bwtraceboundary
    startPt = round(startPt);
    if lineImage(startPt(2), startPt(1)) == 0
        error('Start point must lie on a non-zero pixel in lineImage.');
    end
    
    % Trace boundary from the start point
    boundary = bwtraceboundary(lineImage, [startPt(2), startPt(1)], 'N'); % returns [row, col]
    if isempty(boundary)
        maskSequence = [];
        return;
    end

    % Check which masks the contour passes through
    N = size(masks, 3);
    for m=1:N
        masks(:,:,m)=imdilate(masks(:,:,m)>0, strel('disk', 2));
    end
    maskSequence = zeros(size(boundary, 1), 1);
    for i = 1:size(boundary, 1)
        r = boundary(i, 1);
        c = boundary(i, 2);
        if r >= 1 && r <= size(masks, 1) && c >= 1 && c <= size(masks, 2)
            pixelStack = squeeze(masks(r, c, :));
            idx = find(pixelStack, 1);
            if ~isempty(idx)
                maskSequence(i) = idx;
            end
        end
    end

    % Remove zeros and consecutive duplicates
    maskSequence = maskSequence(maskSequence > 0);
    maskSequence = maskSequence([true; diff(maskSequence) ~= 0]);
    maskSequence = maskSequence(1:ceil(length(maskSequence)/2));
end
