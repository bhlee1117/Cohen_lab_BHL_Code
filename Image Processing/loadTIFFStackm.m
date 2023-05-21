function out = loadTIFFStackm(path, x, y, n)
% function out = loadTIFFStackm(path, x, y, n)
% creates a 3-d array of dimension x, y, n/2, containing the n/2 even numbered
% images in a tiff stack of length n (n must be even).
% Each image has dimension x by y.  We assume that the even numbered images
% are taken in "minus" CPL.
% Yiqiao Tang 7/1/08

out = zeros(x, y, n/2);
for k = 2:2:n;
       out(:,:,k/2) = imread(path, 'tif', k);
       k
end;