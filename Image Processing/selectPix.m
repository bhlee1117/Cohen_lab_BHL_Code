function [Mask, ROI] = selectPix(Img1, Img2, dispImg);
% function [pixMask, ROI] = selectPix(Img1, Img2, dispImg);
%
% Lets a user select one or more regions on a 2D plot of two pixel parameters, and
% produces masks showing the selected pixels
% Img1 and Img2 contain the parameters that will be displayed on the x- and
% y-axes in a scatter plot.  dispImg is the image that will be shown.
%
% ROI is a cell array containing X- and Y-coordinates for each selected ROI.
% Mask is a cell array containing the mask for each selected ROI.
%
% Use with applySelectPix to apply a set of ROIs to a new image pair.
%
% AEC 25 May 2016

figure
clf
colorimg = repmat(mat2gray(dispImg), [1 1 3]);
subplot(1,2,2)
imshow(colorimg, []);

subplot(1,2,1)
plot(Img1(:), Img2(:), '.')
grid on 
grid minor
title('Click the region you want')

k = 1;
Xv = [1 1 1];
while 1
    subplot(1,2,1)
    [Xv, Yv] = getline(gca);
    if length(Xv) < 3;
        break
    end
    hold all
    plot([Xv; Xv(1)], [Yv; Yv(1)]);
    Mask1 = inpolygon(Img1, Img2, Xv, Yv);
    colorDir = rand(1,3);
    colorDir = colorDir/sqrt(sum(colorDir.^2));

    for j = 1:3;
        colorimg(:,:,j) = colorimg(:,:,j) + colorDir(j)*Mask1;
    end;
    subplot(1,2,2)
    imshow(colorimg, [])
    
    ROI{k} = [Xv Yv];
    Mask{k} = Mask1;
    k = k + 1;
end;
hold off



