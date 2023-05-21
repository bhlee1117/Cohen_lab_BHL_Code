function out = loadTIFFStackp2(path, x, y, n)
% function out = loadTIFFStackp2(path, x, y, n)
% creates a 3-d array of dimension x, y, n/2,and here n = 60 as I used in
% image acquisition. The 11-20, 31-40, and 51-60 of the Tiff stack are
% taken under "plus" CPL.
% Yiqiao Tang 4/25/09

out = zeros(x, y, n/2);
for k = 1:10;
       out(:,:,k) = imread(path, 'tif', k+10);
       k
end;
for k = 11:20;
    out(:,:,k) = imread(path,'tif',k+20);
    k
end;
for k=21:30;
    out(:,:,k) = imread(path,'tif',k+30);
    k
end;