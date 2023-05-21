function Mask = applySelectPix(Img1, Img2, dispImg, ROI);
% function Mask = applySelectPix(Img1, Img2, dispImg, ROI);
%
% Apply a set of constraints defined in ROI to the 2-D scatterplot of Img2
% vs. Img1.  Create masks showing the pixels that satisfy the constraints
% and display on dispImg.
% 
% ROI is a cell array containing X- and Y-coordinates for each selected ROI.
% Mask is a cell array containing the mask for each selected ROI.
%
% Use with selectPix to manually define the ROIs.
%
% AEC 25 May 2016

figure
clf
colorimg = repmat(mat2gray(dispImg), [1 1 3]);
subplot(1,2,2)
imshow(colorimg, []);

subplot(1,2,1)
plot(Img1(:), Img2(:), '.')
hold all
nROI = length(ROI);

for j = 1:nROI;
    subplot(1,2,1)
    Xv = ROI{j}(:,1);
    Yv = ROI{j}(:,2);
    plot([Xv; Xv(1)], [Yv; Yv(1)]);
    Mask1 = inpolygon(Img1, Img2, Xv, Yv);
 
    colorDir = rand(1,3);
    colorDir = colorDir/sqrt(sum(colorDir.^2));
    for k = 1:3;
        colorimg(:,:,k) = colorimg(:,:,k) + colorDir(k)*Mask1;
    end;
    subplot(1,2,2)
    imshow(colorimg, [])
    Mask{j} = Mask1;
end;
hold off



