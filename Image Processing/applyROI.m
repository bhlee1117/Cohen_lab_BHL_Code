function [out, bBox] = applyROI(mov, dispImg, ROI);
% [out, bBox] = applyROI(mov, dispImg, ROI);
% returns a sub-region of a movie defined by an input ROI
% mov: input movie.  Can be either a 3D array (Y,X,T) or a 4D color array (Y,X,C,T);
% dispImg: image to display to show ROI location.  If not included, then
% dispImg is the average of the movie along its last dimension.
% ROI: (x,y) coordinates of the ROI
% out: rectangular movie set by the ROI.  Regions outside the
% indicated lines are set to the 20th percentile value within the ROI (over all frames).
% bBox: Rectangular bounding box for the selected region, [xMin yMin Width Height]

nDims = length(size(mov));

figure(99); clf
imshow2(dispImg, []);
hold all
plot([ROI(end, 1); ROI(:,1)], [ROI(end,2); ROI(:,2)]);
hold off

xMin = round(min(ROI(:,1))); xMax = round(max(ROI(:,1)));
yMin = round(min(ROI(:,2))); yMax = round(max(ROI(:,2)));
bBox = [xMin yMin (xMax - xMin) (yMax - yMin)];
if nDims == 3;
    out = mov(yMin:yMax, xMin:xMax,:);
elseif nDims == 4;
    out = mov(yMin:yMax, xMin:xMax,:,:);
end;
out = double(out);
roi2 = [ROI(:,1) - xMin + 1, ROI(:,2) - yMin + 1];
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
        

