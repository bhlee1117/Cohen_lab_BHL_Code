function [images,meanintensities,firstframe] = imagecyclestack_align(fname, imagespercycle, numimages)
% function [images,meanintensities] = imagecyclestack_align(fstem, imagespercycle, numimages)
% Imports a particular number of images from a set of tiff files, averages
% them cyclically, and returns a 3-D array containing the averaged images
% Attempts to align the images in the cycle

img = double(imread(fname, 'tif', 1)); %get the first image to set dimensions
tempimages = zeros([size(img),imagespercycle]);
tempmean = zeros(1,imagespercycle);
for i = 1:imagespercycle %Load the first cycle's images to find alignment point
    fprintf('Image %d\n',i);
    tempimages(:,:,i) = double(imread(fname,'tif',i+1));
    tempmean(i) = mean(mean(tempimages(:,:,i)));
end
avgint = mean(tempmean);
[maxint,maxindex] = max(tempmean);
[minint,minindex] = min(tempmean);
if(maxint - avgint > avgint - minint) %The maximum images deviate most
    threshold = (maxint + avgint)/2;
    firstindex = maxindex;
    if firstindex == 1
        newfirstindex = imagespercycle;
    else
        newfirstindex = firstindex - 1;
    end
    while tempmean(newfirstindex) > threshold
        firstindex = newfirstindex;
        if newfirstindex == 1
            newfirstindex = imagespercycle;
        else
            newfirstindex = newfirstindex - 1;
        end
    end
else %The minimum images deviate most
    threshold = (minint + avgint)/2;
    firstindex = minindex;
    if firstindex == 1
        newfirstindex = imagespercycle;
    else
        newfirstindex = firstindex - 1;
    end
    while tempmean(newfirstindex) < threshold
        firstindex = newfirstindex;
        if newfirstindex == 1
            newfirstindex = imagespercycle;
        else
            newfirstindex = newfirstindex - 1;
        end
    end
end

numcycles = floor((numimages - firstindex)/imagespercycle);

lastindex = firstindex + numcycles * imagespercycle;

images = zeros([size(img),imagespercycle]);
meanintensities = zeros(1,lastindex + 1 - firstindex);
%images(:,:,1) = img;
for i = firstindex:imagespercycle %These images were already loaded - don't reload
    images(:,:,i + 1 - firstindex) = tempimages(:,:,i);
    meanintensities(i) = tempmean(i);
end

for i = imagespercycle+1:lastindex
    fprintf('Image %d\n',i);
    img(:,:) = double(imread(fname,'tif',i+1));
    index = mod(i-firstindex,imagespercycle)+1;
    images(:,:,index) = images(:,:,index) + img(:,:);
    meanintensities(i) = mean(mean(img));
end
for i = 1:imagespercycle
    images(:,:,i) = images(:,:,i) / numcycles;
end

firstframe = firstindex - 1; %Actual index of the first image used from the tif