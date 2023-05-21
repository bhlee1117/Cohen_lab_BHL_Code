function out = radavg(img, x0, y0, r)
% out = radavg(img, x0, y0, r)
% Calculates the radial average of the 2-d matrix img.
% x0 and y0 are the center of the circular average, and 
% need not be integers.  Uses cubic interpolation.
% r is a vector of radii.
% AEC 4 Sept 06

xi = x0 + r'*cos((0:.01:.99)*2*pi);
yi = y0 + r'*sin((0:.01:.99)*2*pi);
zi = interp2(img, xi, yi, 'cubic')';
out = mean(zi);
