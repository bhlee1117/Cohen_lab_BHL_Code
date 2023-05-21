function out = loadTIFFs(path, x, y, n)
% function out = loadTIFFs(path, x, y, n)
% creates a 3-d array of dimension x, y, n, containing the n images in a
% tiff stack, where each image has dimension x by y.
% AEC 9/15/05
% JHH 25 March 09 Suppressed frame text printout
% SL 6 April 09 Output as double

out = zeros(x, y, n);
for k = 1:n;
    out(:,:,k) = double(imread(path, 'tif', k));
end;