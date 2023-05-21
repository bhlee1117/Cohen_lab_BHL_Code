function out = binTIFFStack(path, bin, n)
% function out = binTIFFStack(path, bin, n)
% creates a 3-d array containing the n images in a
% tiff stack, where each image has been subjected to bin x bin binning.
% AEC 4/28/09


a = imfinfo(path);
xsize = a(1).Width;
ysize = a(1).Height;
out = zeros(ysize/bin, xsize/bin, n);
x2 = length(1:bin:xsize);
y2 = length(1:bin:ysize);
for k = 1:n;
    k
    tmp = double(imread(path, 'tif', k));
    tmp2 = zeros(y2, x2);
    for r = 1:bin;
        for c = 1:bin;
            tmp2 = tmp2 + tmp(r:bin:end-bin+r, c:bin:end-bin+c);
        end;
    end;
    out(:,:,k) = tmp2/(bin^2);
end;