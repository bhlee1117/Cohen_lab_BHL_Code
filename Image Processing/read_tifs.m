function out = read_tifs(filename, firstimage, lastimage)
% function out = read_tifs(filename, firstimage, lastimage)
% Reads a series of tiff images into a single 3-d array
% AEC 4/17/05


if firstimage < 10;
    img1 = imread([filename '000' num2str(firstimage) '.tif'], 'tif');        
elseif firstimage < 100;
    img1 = imread([filename '00' num2str(firstimage) '.tif'], 'tif');     
elseif firstimage < 1000;
    img1 = imread([filename '0' num2str(firstimage) '.tif'], 'tif');     
else
    img1 = imread([filename num2str(firstimage) '.tif'], 'tif');     
end;

[rows, cols] = size(img1);
out = uint16(zeros(lastimage - firstimage + 1, rows, cols));

for j = firstimage:lastimage
    if j < 10;
        out(j-firstimage + 1,:,:) = imread([filename '000' num2str(j) '.tif']);        
    elseif j < 100;
        out(j-firstimage + 1,:,:) = imread([filename '00' num2str(j) '.tif']);     
    elseif j < 1000;
        out(j-firstimage + 1,:,:) = imread([filename '0' num2str(j) '.tif']);     
    else
        out(j-firstimage + 1,:,:) = imread([filename num2str(j) '.tif']);     
    end;
end;
