%Load all images from .sif file into a 512x512x50 matrix
function out = andortest(filename);
for n=1:50, stack(:,:,n)=readCCDim(filename,2,n); end

%Calculate average intensity by averaging a single pixel through time
averageintensity = mean(stack,3);
%View image
imview(averageintensity)

%Calculate variance of a pixel through time
varianceofpixels = var(stack,0,3);
%View image
imview(varianceofpixels)

%Reshape matrix into a big vector
reshapedstack = reshape(stack,13107200,1);

%Calculate histogram (needs tweaking because x isn't really working)
x=0:1:500;
hist(reshapedstack,x)