function mask=refineMask(BW,tmpImg)

mask=BW;
mask=imfill(bwareaopen(BW,30,'holes');
tmpImg=imhimin(rgb2gray(tmpImg),13);
tmpImg=watershed(tmpImg);
tmpImg=tmpImg==0;
tmpImg=bwareaopen(tmpImg,200,8);
mask(tmpImg)=0;