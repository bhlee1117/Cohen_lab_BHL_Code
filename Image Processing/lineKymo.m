function kymo = lineKymo(mov, dL, dP, displayImg);
% function kymo = lineKymo(mov, dL, dP, displayImg);
% Takes a movie 'mov' and displays the image 'displayImg'.
% If 'displayImg' is not specified, displays the mean of the movie along
% the third axis.
% User inputs one line.  The line is divided into rectangular segments of length
% 'dL' and perpendicular width 'dP'.  The average of each rectangle is given as a function
% of time.
% AEC 8 March 2021

% Check the number of arguments
if nargin == 3;
    displayImg = mean(mov, 3);
end;
figure(1); clf
imshow(displayImg, [])

% Get user to specify the line to average along
title('Left click to start line; right click to end')
while 1
    [X, Y] = getline(gcf);
    if length(X) > 2
        title('Error, line can only have two points')
    else
        break
    end;
end;
hold all
plot(X, Y, 'r-')
hold off

dX = X(2) - X(1); % Vector pointing from start to end of line
dY = Y(2) - Y(1);

dX0 = dX/(dX^2 + dY^2)^.5; % Unit vector along line
dY0 = dY/(dX^2 + dY^2)^.5;

dXp = -dY0; % Unit vector perpendicular to line
dYp = dX0;

L = (dX^2 + dY^2)^.5;  % Length of the line

lSteps = 0:dL:L;  % displacements along the line
nSteps = length(lSteps);  % total number of steps

% Define the coordinates of all the rectangles
clear roi
for j = 1:nSteps;
    Xc = X(1) + lSteps(j)*dX0; % box center
    Yc = Y(1) + lSteps(j)*dY0;
    Xs = [Xc - dX0*dL/2 - dXp*dP/2; % Make closed paths (5 points)
          Xc - dX0*dL/2 + dXp*dP/2; 
          Xc + dX0*dL/2 + dXp*dP/2;
          Xc + dX0*dL/2 - dXp*dP/2;
          Xc - dX0*dL/2 - dXp*dP/2];
    Ys = [Yc - dY0*dL/2 - dYp*dP/2; 
          Yc - dY0*dL/2 + dYp*dP/2; 
          Yc + dY0*dL/2 + dYp*dP/2; 
          Yc + dY0*dL/2 - dYp*dP/2;
          Yc - dY0*dL/2 - dYp*dP/2];
      roi{j} = [Xs, Ys];
end;
kymo = apply_clicky(roi, mov, 'no');
