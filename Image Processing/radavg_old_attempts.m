function out = radavg(img, x0, y0, r)
% out = radavg(img, x0, y0)
% Calculates the radial average of the 2-d matrix img.
% x0 and y0 are the center of the circular average, and 
% need not be integers.
% Works better if size(img) is odd.
%THIS ROUTINE MIGHT BE BROKEN AEC 1/6/06

% The calculation weights the amplitude in each grid-element
% by the arc-length of the circle going through that grid-element
% AEC 29 April 2005


if 0;
m = size(img, 1);
if round(m/2) == m/2;
    out = zeros((m/2) + 1,1);
else;
    out = zeros((m+1)/2,1);
end;

T = maketform('affine', [1 0 0; 0 1 0; 0 0 1]);
x = x0 - (m+1)/2;
y = y0 - (m+1)/2;

imgctr = imtransform(img, T, 'bicubic', 'XData', [1+x,m+x], 'YData', [1+y,m+y]);
c = 0;
for j = 0:5:360;
    c = c + 1;
    tmp = imrotate(imgctr, j, 'bicubic');
    if round(m/2) == m/2;
        tmp = imtransform(img, T, 'bicubic', 'XData', [1.5,m+.5], 'YData', [1.5,m+.5]);
        out = out + tmp(m/2, m/2:m)';
    else;
        out = out + tmp((m+1)/2,((m+1)/2):m)';
    end;
end;
out = out/c;
end; % if 0


if 0;
pcolor(img);    % plot the image
shading 'interp'
hold on
nr = max(size(r));
out = zeros(nr,1);

xs = 0;
ys = 0;

for j = 1:nr;
    if (ceil(y0 - r(j)) <= floor(y0 + r(j)));  % the circle crosses at least one line y = integer
        yints = ceil(y0-r(j)):floor(y0+r(j));   % values of y where the circle crosses the lines y = integer
        xplus = x0 + sqrt(r(j)^2 - (yints - y0).^2);  % x-coordinates of the intersections
        xminus = x0 - sqrt(r(j)^2 - (yints - y0).^2);
        ys = max(size(yints));
    end;
    if (ceil(x0 - r(j)) <= floor(x0 + r(j)));  % the circle crosses at least one line x = integer
        xints = ceil(x0-r(j)):floor(x0+r(j));   % values of x where the circle crosses the lines x = constant
        yplus = y0 + sqrt(r(j)^2 - (xints - x0).^2);   % y-coordinates of the intersections
        yminus = y0 - sqrt(r(j)^2 - (xints - x0).^2);
        xs = max(size(xints));
    end;
    
    if xs && ys;  % the circle intersects x and y lines
        intcepts = [xplus', yints'; xminus', yints'; xints', yplus'; xints', yminus'];  % an unordered list of the intersections
    elseif ys;  % only intersects y-lines
        intcepts = [xplus', yints'; xminus', yints'];
    elseif xs;  % only intersects x-lines
        intcepts = [xints', yplus'; xints', yminus']; 
    end;
    
    if xs | ys;  % the circle has >= 1 intersection
        numints = max(size(intcepts));  % count the number of intersections
        thetas = atan2(intcepts(:,2) - y0, intcepts(:,1)- x0);  % find the angle at each intersection
        [thetas, ind] = sort(thetas); % sort the angles ascending from -pi to pi
        
        % this part is only needed for the plotting function
        for k = 1:numints
            intcepts2(k,:) = intcepts(ind(k),:);  % apply the sorting to the x-y coordinates of each intersection
        end;
        % ---------------------------------------------
        
        for k = 1:(numints-1)   % weight each contribution from img by the arc-length
            xmidpt = x0 + r(j)*cos((thetas(k) + thetas(k+1))/2);
            ymidpt = y0 + r(j)*sin((thetas(k) + thetas(k+1))/2);
            out(j) = out(j) + img(floor(ymidpt), floor(xmidpt))*(thetas(k+1) - thetas(k));
        end;
        % the last arc connecting the beginning and end is handled
        % separately
        xmidpt = x0 + r(j)*cos((thetas(1) + 2*pi + thetas(numints))/2);
        ymidpt = y0 + r(j)*sin((thetas(1) + 2*pi + thetas(numints))/2);
        out(j) = out(j) + img(floor(ymidpt), floor(xmidpt))*(thetas(1) + 2*pi - thetas(numints));
        out(j) = out(j)/(2*pi);
        
        plot(intcepts2(:,1), intcepts2(:,2))
    else
        out(j) = img(floor(y0), floor(x0));
    end;
end;
hold off
end; %if 0;