function out = loadTIFFStack(path, x, y, n)
% function out = loadTIFFStack(path, y, x, n)
% creates a 3-d array of dimension y, x, n, containing the n images in a
% tiff stack, where each image has dimension x by y.
% AEC 9/15/05

out = zeros(x, y, n);
for k = 1:n;
    if ~mod(k, 10);
        k
    end;
    out(:,:,k) = imread(path, 'tif', k);
end;