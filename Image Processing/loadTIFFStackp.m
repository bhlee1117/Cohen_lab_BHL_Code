function out = loadTIFFStackp(path, x, y, n)
% function out = loadTIFFStackp(path, x, y, n)
% creates a 3-d array of dimension x, y, n/2, containing the n/2 odd numbered
% images in a tiff stack of length n (n must be even).
% Each image has dimension x by y.  We assume that the odd numbered images
% are taken in "plus" CPL.
% Yiqiao Tang 7/1/08

out = zeros(x, y, n/2);
for k = 1:2:(n-1);
       out(:,:,(k+1)/2) = imread(path, 'tif', k);
       k
end;