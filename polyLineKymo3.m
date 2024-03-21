function [kymo,kymoROI] = polyLineKymo3(mov, dL, dP, displayImg);
% function [kymo, kymoROI] = polyLineKymo2(mov, dL, dP, displayImg);
% Takes a movie 'mov' and displays the image 'displayImg'.
% If 'displayImg' is not specified, displays the mean of the movie along
% the third axis.
% User inputs a line of arbitrary number of segments.  The line is divided into rectangular segments of length
% 'dL' and perpendicular width 'dP'.  The average of each rectangle is given as a function
% of time.
% kymoROI: the X, Y coordinates of the clicked spots.
% modified from lineKymo.m which only takes straight lines
% modified from polyLineKymo.m which doesn't return the clicked
% coordinates.
% AEC 8 March 2021

% Check the number of arguments
if nargin == 3;
    displayImg = mean(mov, 3);
end;
figure(1); clf
imshow(displayImg, [])
[ySize, xSize] = size(displayImg);

% Get user to specify the line to average along
title('Left click to start line; right click to end')
[X, Y] = getline(gcf);

hold all
plot(X, Y, 'r-')
hold off

dX = diff(X); % Vector pointing from start to end of line
dY = diff(Y);
nSeg = length(dX); % number of line segments

dX0 = dX./(dX.^2 + dY.^2).^.5; % Unit vector along line
dY0 = dY./(dX.^2 + dY.^2).^.5;

dXp = -dY0; % Unit vector perpendicular to line
dYp = dX0;

L = (dX.^2 + dY.^2).^.5;  % Length of the line
c = 1;
clear roi
trail=[];
for k = 1:nSeg;
    lSteps = 0:dL:L(k);  % displacements along the line
    nSteps = length(lSteps);  % total number of steps

    % Define the coordinates of all the rectangles
    for j = 1:nSteps;
        Xc = X(k) + lSteps(j)*dX0(k); % box center
        Yc = Y(k) + lSteps(j)*dY0(k);
        trail=[trail; [Xc Yc]];
        Xs = [Xc - dX0(k)*dL/2 - dXp(k)*dP/2; % Make closed paths (5 points)
              Xc - dX0(k)*dL/2 + dXp(k)*dP/2; 
              Xc + dX0(k)*dL/2 + dXp(k)*dP/2;
              Xc + dX0(k)*dL/2 - dXp(k)*dP/2;
              Xc - dX0(k)*dL/2 - dXp(k)*dP/2];
        Ys = [Yc - dY0(k)*dL/2 - dYp(k)*dP/2; 
              Yc - dY0(k)*dL/2 + dYp(k)*dP/2; 
              Yc + dY0(k)*dL/2 + dYp(k)*dP/2; 
              Yc + dY0(k)*dL/2 - dYp(k)*dP/2;
              Yc - dY0(k)*dL/2 - dYp(k)*dP/2];
          Xs(Xs > xSize) = xSize; Xs(Xs < 1) = 1;  % Make sure the rectangles don't go off the image.
          Ys(Ys > ySize) = ySize; Ys(Ys < 1) = 1;
          roi{c} = [Xs, Ys];
          c = c + 1;
    end;
end;
kymo = apply_clicky(roi, mov);  % add a 'no' as the third argument to not show the results.
kymoROI = trail;
