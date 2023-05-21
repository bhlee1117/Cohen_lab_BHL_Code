function [h, outImgs] = getFibers(inImgs, dW, n)
% [h, outImgs] = getFibers(inImgs, dW, n)
% Delineate the fibers in a picture with many fibers.
% inImgs: multicolor fiber image
% dW: width to ascribe to each fiber
% n: Figure number to display the whole image
% h: handle to image with multiple ROI lines
% outImgs: cell array with images of the individual fibers
% Need to draw one more line after hitting a key to terminate.
%
% AEC 4 Sept 2021

[ySize, xSize, nCol] = size(inImgs);
if nCol ~= 3
    disp('Error: need 3 color channels')
    return
end;

figure(n); clf;
imshow2(inImgs)
title('Draw ROIs, press any key to finish')
KEY_IS_PRESSED = 0;
set(gcf, 'KeyPressFcn', @myKeyPressFcn)

c = 1;
while ~KEY_IS_PRESSED
    f = drawline(gca, 'LineWidth', 1);
    if ~KEY_IS_PRESSED
        text(f.Position(1,1), f.Position(1,2), num2str(c), 'Color', 'w', 'FontSize', 16)
        c = c+1;
    end;
end;
delete(f)% get rid of the last line element, which is empty

h = gca;
roi = [];
for j = 1:length(h.Children)  % find all the line ROIs
    if strcmp(h.Children(j).Type, 'images.roi.line')
        roi = [roi h.Children(j)];
    end;
end;
nROI = length(roi);

for j = nROI:-1:1
    if(isempty(roi(j).Position))
        roi(j) = [];
    end;
end;
nROI = length(roi);

figure(33); clf
for j = 1:nROI
    coords = roi(j).Position;  % coordinates of the line
    X = coords(:,1);
    Y = coords(:,2);
    % geometry to rotate the image
    dX = diff(X);
    dY = diff(Y);
    theta = atan2d(dY, dX);
    img2 = imrotate(inImgs, theta);
    [ySize2, xSize2, ~] = size(img2);
    
    % geometry to rotate the line
    X0 = xSize/2; Y0 = ySize/2;
    xRel = X - X0;
    yRel = Y - Y0;
    rotMat = [cosd(-theta) -sind(-theta); sind(-theta) cosd(-theta)];
    coordsF = rotMat*([xRel'; yRel']);
%     XF = coordsF(1,:) + X0;
%     YF = coordsF(2,:) + Y0;
    XF = coordsF(1,:) + xSize2/2;
    YF = coordsF(2,:) + ySize2/2;
    XF = max(XF, 1); YF = max(YF, 1);
    XF = min(XF, xSize2); YF = min(YF, ySize2);

    yBot = max(YF(1)-dW, 1);
    yTop = min(YF(1)+dW, ySize2);
    outImgs{j} = img2(round(yBot:yTop),round(XF(1):XF(2)),:);
    for k = 1:3
        outImgs{j}(:,:,k) = mat2gray(outImgs{j}(:,:,k));
    end;
    subplot(nROI,1,j)
    imshow2(outImgs{j});
end;


function myKeyPressFcn(hObject, event)
KEY_IS_PRESSED  = 1;
end


end
