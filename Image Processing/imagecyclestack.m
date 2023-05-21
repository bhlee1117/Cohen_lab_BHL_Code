function [images,meanintensities] = imagecyclestack(fname, imagespercycle, numimages)
% function images = imagecycle(fstem, imagespercycle, numimages)
% Imports a particular number of images from a set of tiff files, averages
% them cyclically, and returns a 3-D array containing the averaged images

img = double(imread(fname, 'tif', 1)); %get the first image to set dimensions
images = zeros([size(img),imagespercycle]);
meanintensities = zeros(1,numimages-1);
%images(:,:,1) = img;
numimgs = zeros(imagespercycle);
for i = 1:numimages-1
    fprintf('Image %d\n',i);
    img(:,:) = double(imread(fname,'tif',i+1));
    images(:,:,mod(i-1,imagespercycle)+1) = images(:,:,mod(i-1,imagespercycle)+1) + img(:,:);
    numimgs(mod(i-1,imagespercycle)+1) = numimgs(mod(i-1,imagespercycle)+1) + 1;
    meanintensities(i) = mean(mean(img));
end
for i = 1:imagespercycle
    images(:,:,i) = images(:,:,i) / numimgs(i);
end