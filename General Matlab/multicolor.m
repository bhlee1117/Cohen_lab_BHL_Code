function h = multicolor(imgs, colorMat);
% function h = multicolor(imgs, colorMat);
% display a multicolor image of an image stack with an arbitrary number of
% frames.
% Each frame is automatically scaled min:max.
% colorMat is an optional nx3 matrix (n = number of imgs) specifying the
% color value of each image.  If omitted, the colors are taken from the
% 'jet' colormap.
%
% Adam Cohen 16 March 2022

[ySize, xSize, nCol] = size(imgs);
if nargin == 1
    colorMat = generateColorSpecLocal(nCol);
end;

colorImg = zeros(ySize, xSize, 3);
for j = 1:nCol;
    for k = 1:3;
        colorImg(:,:,k) = colorImg(:,:,k) + colorMat(j,k)*mat2gray(imgs(:,:,j));
    end;
end
imshow2(colorImg);
h = gcf;
    