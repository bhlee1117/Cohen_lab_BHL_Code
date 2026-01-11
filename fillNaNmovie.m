function strImgbkg_filled=fillNaNmovie(strImgbkg)

strImgbkg_filled=[];
for n=1:size(strImgbkg,3)
    A=strImgbkg(:,:,n);
[~, idx]=bwdist(~isnan(A));
A(isnan(A)) = A(idx(isnan(A)));  % assign nearest values
strImgbkg_filled(:,:,n)=A;
end