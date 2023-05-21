function out = read_tifs3(filename, firstimage, lastimage)
% function out = read_tifs3(filename, firstimage, lastimage)
% Reads a series of tiff images into a single 3-d array
% AEC 4/17/05
% Modified 7/7/08 AEC and Malcolm Campbell
% Assumes filenames do not have leading zeros before the image number.


img1 = imread([filename num2str(firstimage) '.tif'], 'tif');     

[rows, cols] = size(img1);
out = uint16(zeros(rows, cols, lastimage - firstimage + 1));

for j = firstimage:lastimage
    out(:,:,j-firstimage + 1) = imread([filename num2str(j) '.tif'], 'tif');     
end;
out = double(out);
