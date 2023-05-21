function kymo = wedgeKymo(in, nSegments, displayImg);
% function kymo = wedgeKymo(in, nSegments);
% Takes a movie 'in' and displays the image 'displayImg'.
% If 'displayImg' is not specified, displays the mean of the movie along
% the third axis.
% User inputs two lines.  The lines are divided into segments of length
% 'nSegments', and the average of each quadrilateral is given as a function
% of time.
% AEC 24 May 2016

if nargin == 2;
    displayImg = mean(in, 3);
end;
figure(1); clf
subplot(1,2,1);
imshow(displayImg, [])

title('Click a line for the first boundary')
while 1
    [x1, y1] = getline(gcf);
    if length(x1) > 2
        title('Error, line can only have two points')
    else
        break
    end;
end;
hold all
plot(x1, y1, 'r-')

title('Click a line for the second boundary')
while 1
    [x2, y2] = getline(gcf);
    if length(x2) > 2
        title('Error, line can only have two points')
    else
        break
    end;
end;
hold all
plot(x2, y2, 'r-')

x1Vec = linspace(x1(1), x1(2), nSegments + 1);
y1Vec = linspace(y1(1), y1(2), nSegments + 1);
x2Vec = linspace(x2(1), x2(2), nSegments + 1);
y2Vec = linspace(y2(1), y2(2), nSegments + 1);

plot([x1Vec; x2Vec], [y1Vec; y2Vec], 'c-');

[ySize, xSize, nFrames] = size(in);
kymo = zeros(nSegments, nFrames);
XCoords = repmat(1:xSize, [ySize, 1]);
YCoords = repmat((1:ySize)', [1, xSize]);
for j = 1:nSegments;
    inMat = find(inpolygon(XCoords, YCoords, [x1Vec(j), x1Vec(j+1), x2Vec(j+1), x2Vec(j)], [y1Vec(j), y1Vec(j+1), y2Vec(j+1), y2Vec(j)]));
    for k = 1:nFrames;
        pixIdx = inMat + (k-1)*ySize*xSize;
        kymo(j,k) = mean(in(pixIdx));
    end;
    j
end;
subplot(1,2,2);
pcolor(kymo'); shading 'interp'


