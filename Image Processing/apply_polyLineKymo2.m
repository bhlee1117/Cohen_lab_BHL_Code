function kymo = apply_polyLineKymo2(mov, dL, dP, displayImg, kymoROI);
% function kymo = apply_polyLineKymo2(mov, dL, dP, displayImg, kymoROI);
% Takes a movie 'mov' and displays the image 'displayImg'.
% If 'displayImg' is not specified, displays the mean of the movie along
% the third axis.
% kymoROI is an [X, Y] line of arbitrary number of segments.  The line is divided into rectangular segments of length
% 'dL' and perpendicular width 'dP'.  The average of each rectangle is given as a function
% of time.

% Check the number of arguments
if nargin == 3;
    displayImg = mean(mov, 3);
end;
figure(1); clf
imshow(displayImg, [])
[ySize, xSize] = size(displayImg);

% Coordinates of the line
X = kymoROI(:,1);
Y = kymoROI(:,2);

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
for k = 1:nSeg;
    lSteps = 0:dL:L(k);  % displacements along the line
    nSteps = length(lSteps);  % total number of steps

    % Define the coordinates of all the rectangles
    for j = 1:nSteps;
        Xc = X(k) + lSteps(j)*dX0(k); % box center
        Yc = Y(k) + lSteps(j)*dY0(k);
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
          Xs(Xs > xSize) = xSize; Xs(Xs < 1) = 1;
          Ys(Ys > ySize) = ySize; Ys(Ys < 1) = 1;
          roi{c} = [Xs, Ys];
          c = c + 1;
    end;

end;
kymo = apply_clicky(roi, mov);

