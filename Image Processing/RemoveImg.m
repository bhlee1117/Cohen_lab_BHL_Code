function [out, scaleVecs] = RemoveImg(mov, refImg, remOffset); 

% function [out, scaleVecs] = RemoveImg(mov, refImg, remOffset); 
% project out one or more images from a movie.
% size(Movie) = [Y, X, T];
% size(refImg) = [Y, X, N], where N is the number of reference images to
% remove.
% if refImg == 0, then only the average intensity is subtracted from each
% frame.
% 
% remOffset: Boolean variable.  If 1, add an images of 1's to the refImg, to
% remove constant uniform background from each frame.  If 0, do not.  If omitted,
% default = 1.
%
% out = corrected movie
% scaleVecs = vectors of the amplitude coefficients for each image removed.
% size(scaleVecs) = [T, N], or [T, N+1] if remOffset == 1.
%
% AEC + LF 5 July 2016


if nargin == 2;
    remOffset = 1;  % Default: if not specified whether to remove offset, then subtract off mean intensity.
end;
if ~isscalar(remOffset)
    'remOffset must be a scalar (zero or nonzero)'
    return
end;
% 
[ySize, xSize, nFrame] = size(mov);

if refImg == 0;  % Then only subtract the average intensity from each frame.
    intens = mean(mean(mov));
    out = mov - repmat(intens, [ySize, xSize, 1]);
    scaleVecs = intens;
else;  %Project out the images in refImg.
    % Make sure refImg is the right size and orientation
    [refY, refX, refN] = size(refImg);
    if ~(refY == ySize & refX == xSize);
        'Size of refImg must match x- and y- size of movie'
        return
    end;
    
    if remOffset;
        dRef = refImg - repmat(mean(mean(refImg)), [ySize, xSize, 1]); % remove the DC offset from refSig because we will account for offset separately.
        % Add a row of ones to measure the DC offset at each pixel.
        I = zeros(ySize, xSize, refN + 1);
        I(:,:,1) = 1;
        I(:,:,2:end) = dRef;
    else
        I = refSig;
    end;

    IVec = tovec(I);
    movVec = tovec(mov);
   
    C = inv(IVec'*IVec)*IVec'*movVec;  % Linear algebra to find the pseudo-inverse
    scaleVecs = C';  % convert weights to images.

    residMovVec = movVec - IVec*C;  % Look at the residuals
    out = toimg(residMovVec, ySize, xSize);
end;


