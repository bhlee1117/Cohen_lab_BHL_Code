function out = loadTIFFIntens(path, n)
% function out = loadTIFFIntens(path, n)
% Read n frames of a Tiff stack and return 
% the intensity of each frame.

out = zeros(1, n);
for k = 1:n;
    k
    tmp = imread(path, 'tif', k);
    out(k) = mean(mean(tmp));
end;