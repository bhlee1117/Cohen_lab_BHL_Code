function mov = outprod(img, timeseries);

% function mov = outprod(img, timeseries);
% Takes a 2-D image and produces a 3-D array:
% mov(:,:,j) = img*timeseries(j);
%
% AEC 23 May 2015

[ySize, xSize] = size(img);
nFrames = length(timeseries);
mov = zeros(ySize, xSize, nFrames);
for j = 1:nFrames;
    mov(:,:,j) = img.*timeseries(j);
end;
