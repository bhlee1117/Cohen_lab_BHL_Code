function [meanval, stddev, npix, pixvals] = clickroi(dataimg, visimg)
% function [meanval, stddev, npix, pixvals] = clickroi(dataimg, visimg);
% Returns statistics on the pixels within a user-selected ROI on an image.
% dataimg is the image to which the calculations are applied
% visimg is the image displayed to the user.
% visimg and dataimg must be the same size.
% If only one image is supplied, then the same image is used for both
% purposes
% AEC 26 Feb. 2012

if nargin == 1
    visimg = dataimg;
end;

if size(dataimg) ~= size(visimg)
    'Error: Images must be the same size'
    return
end;

figure
imshow(visimg, [], 'InitialMagnification', 'fit')
hold on;

[ysize, xsize] = size(dataimg);
[x, y] = meshgrid(1:xsize, 1:ysize);
[xv, yv] = (getline(gca, 'closed'));
inpoly = inpolygon(x,y,xv,yv);

plot(xv, yv, 'r-', 'Linewidth', 1);
hold off

npix = sum(inpoly(:));
pixvals = dataimg(inpoly(:));
meanval = mean(pixvals);
stddev = std(pixvals);
