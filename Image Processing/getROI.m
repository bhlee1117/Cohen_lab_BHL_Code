function [out, bBox, roi] = getROI(mov, dispImg);
% [out, bBox, roi] = getROI(mov, dispImg);
% returns a sub-region of a movie defined by freehand drawing
% mov: input movie.  Can be either a 3D array (Y,X,T) or a 4D color array (Y,X,C,T);
% dispImg: image to display for selecting ROIs.  If not included, then
% dispImg is the average of the movie along its last dimension.
% out: rectangular movie set by the ROI.  Regions outside the
% indicated lines are set to the 20th percentile within the ROI (over all frames).
% roi: (x,y) coordinates of the ROI
% bBox: Rectangular bounding box for the selected region, [xMin yMin Width Height]

nDims = length(size(mov));
if nargin == 1;
    dispImg = mean(dispImg, nDims);
end;

figure(99); clf
imshow2(dispImg, []);
title('Draw the ROI, hit any key when done')
r = drawfreehand(gca);
'Hit any key when done';
pause

roi = r.Position;
xMin = round(min(roi(:,1))); xMax = round(max(roi(:,1)));
yMin = round(min(roi(:,2))); yMax = round(max(roi(:,2)));
bBox = [xMin yMin (xMax - xMin) (yMax - yMin)];
if nDims == 3;
    out = mov(yMin:yMax, xMin:xMax,:);
elseif nDims == 4;
    out = mov(yMin:yMax, xMin:xMax,:,:);
end;
out = double(out);
roi2 = [roi(:,1) - xMin + 1, roi(:,2) - yMin + 1];
[ySize, xSize, nFrames] = size(out);
[X, Y] = meshgrid(1:xSize, 1:ySize);
inside = inpolygon(X, Y, roi2(:,1), roi2(:,2));
if nDims == 3;
    fMin = prctile(out(:), 20);
    out = out.*inside;
    out = out + (1-inside)*fMin;
elseif nDims == 4;
    for j = 1:3;
        tmp = squeeze(out(:,:,j,:));
        fMin = prctile(tmp(:), 20);
        tmp = tmp.*inside;
        tmp = tmp + (1-inside)*fMin;
        out(:,:,j,:) = reshape(tmp, [ySize, xSize, 1, nFrames]);
    end;
end;
        

