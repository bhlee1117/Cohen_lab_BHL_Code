clear

for i=1:15
im(:,:,i)=imread(['D:\' '1-1.tif'],4*i-3);
end

imfilt=double(imgaussfilt(im(:,:,1),5));
imfilt=imfilt*255/max(reshape(imfilt,size(imfilt,1)*size(imfilt,2),1));
imfilt=uint8(imfilt);

imfilt=2^8-1-imfilt;

circleFinder(imfilt)
