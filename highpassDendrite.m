function strImg=highpassDendrite(ref_im,SWC)

sz=size(ref_im);
[xx,yy] = meshgrid(1:sz(1),1:sz(2));
mask = false(sz(2),sz(1));
r=ones(size(SWC,1),1)*3;
r(1)=10;
for i = 1:size(SWC,1)
    mask = mask | ((xx - SWC(i,1)).^2 + (yy - SWC(i,2)).^2 <= r(i)^2);
end

strImg=ref_im-100;%./imgaussfilt(Result.ref_im,50);
strImgMasked=strImg;
strImgMasked(mask==1)=NaN;
strImgbkg = medfilt2nan(strImgMasked, ones(1,2)*35);
[~, idx]=bwdist(~isnan(strImgbkg));
strImgbkg(isnan(strImgbkg)) = strImgbkg(idx(isnan(strImgbkg)));  % assign nearest values
strImgbkg = imgaussfiltnan(strImgbkg, 5);
strImg=strImg-strImgbkg;
end