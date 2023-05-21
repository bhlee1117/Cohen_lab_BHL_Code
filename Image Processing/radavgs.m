function [out, r] = radavgs(img, x0, y0, r)
% out = radavg(img, x0, y0, r)
% Calculates the radial average of the 2-d matrix img.
% x0 and y0 are the center of the circular average, and 
% need not be integers.  Uses cubic interpolation.
% r is a vector of radii.
% AEC 4 Sept 06
% Vicente Parot 2018

if ~exist('y0','var')
    y0 = size(img,1)/2+.5;
end

if ~exist('x0','var')
    x0 = size(img,2)/2+.5;
end

if ~exist('r','var')
    maxr = max(sqrt(sum(([1 1; size(img)] - [y0 x0]).^2,2)));
    r = linspace(-maxr,maxr,ceil(max(size(img))*2*sqrt(2)));
end

xi = x0 + r'*cos((0:.01:.99)*2*pi);
yi = y0 + r'*sin((0:.01:.99)*2*pi);
for it = 1:size(img,3)
    zi = interp2(img(:,:,it), xi, yi)';
    out(it,:) = nanmean(zi);
end
