function out = loadTIFFStackRange(path, x, y, iStart, iEnd)
% function out = loadTIFFStackRange(path, x, y, iStart, iEnd)
% creates a 3-d array of dimension x, y, n, containing the iEnd-iStart+1 images in a
% tiff stack, where each image has dimension x by y.
% AEC 05 Sept 15
% JHH 25 March 09 Suppressed frame text printout
% JHH 18 March 09 Modified to be able to read a substack with arbitrary
% start and end frames.

out = zeros(x, y, iEnd-iStart+1);

for iFrame = iStart:iEnd
    out(:, :, iFrame-iStart+1) = imread(path, 'tif', iFrame);
end