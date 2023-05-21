function images = imagecycle(fstem, imagespercycle, numimages)
% function images = imagecycle(fstem, imagespercycle, numimages)
% Imports a particular number of images from a set of tiff files, averages
% them cyclically, and returns a 3-D array containing the averaged images

img = double(imread([fstem '_0000.tif'], 'tif')); %get the first image to set dimensions
images = zeros([size(img),imagespercycle]);
%images(:,:,1) = img;
numimgs = zeros(imagespercycle);
for i = 1:numimages-1
    fprintf('Image %d\n',i);
    img(:,:) = double(imread(sprintf('%s_%04.0f.tif',fstem,i),'tif'));
    images(:,:,mod(i-1,imagespercycle)+1) = images(:,:,mod(i-1,imagespercycle)+1) + img(:,:);
    numimgs(mod(i-1,imagespercycle)+1) = numimgs(mod(i-1,imagespercycle)+1) + 1;
end
for i = 1:imagespercycle
    images(:,:,i) = images(:,:,i) / numimgs(i);
end