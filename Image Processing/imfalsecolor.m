function cImage = imfalsecolor(bwImage,colorSel,normSelect)

% function cImage = imfalsecolor(bwImage,colorSel)
% IMFALSECOLOR takes an image, bwimage, and converts it to a image of equal
% size with three dimensions such that it can be displayed in color.  The
% image is normalized to intensinties between 0 and 1 so that it can be
% displayed easily in MATLAB.  Current choices for colorSel are 'blue',
% 'red', 'green', and 'orange'.  The default colorSel is blue.  The program
% only nomralizes the data if normSelect is set to TRUE.

if nargin < 2
    colorSel = 'blue';
    normSelect = false;
elseif nargin < 3
    normSelect = false;
end

colorSel = lower(colorSel);

switch colorSel
    case {'blue'}
        cMap = [0 0 1];
    case {'green'}
        cMap = [0 1 0];
    case {'red'}
        cMap = [1 0 0];
    case {'orange'}
        cMap = [1 .5471 .2018];
    otherwise
        disp('Unknown color selection, reverting to blue.');
        cMap = [0 0 1];
end

[ySize,xSize,numImages] = size(bwImage);
cImage = zeros(ySize,xSize,3,numImages);

for k = 1:numImages
    for i = 1:3
        cImage(:,:,i,k) = cMap(i)*bwImage(:,:,k);
    end
end

if normSelect == true
    cImage = cImage - min(cImage(:));
    cImage = cImage./(max(cImage(:)));
end
end