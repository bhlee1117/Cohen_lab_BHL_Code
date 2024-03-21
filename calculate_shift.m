function [yoffSet xoffSet]=calculate_shift(avgImg,avgImg2)
sz2=size(avgImg2);
sz1=size(avgImg);
crop_avgImg2=avgImg2(sz2(1)*1/4:sz2(1)*3/4,sz2(2)*1/4:sz2(2)*3/4);

c = normxcorr2(mat2gray(crop_avgImg2), mat2gray(avgImg));
[ypeak, xpeak] = find(c==max(c(:)));

yoffSet = ypeak-size(crop_avgImg2,1)/2-sz1(1)/2+0.5;
xoffSet = xpeak-size(crop_avgImg2,2)/2-sz1(2)/2+0.5;