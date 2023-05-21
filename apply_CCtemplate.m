function xcorr_img_th=apply_CCtemplate(Result_perm,Patt,corr_th)
wind_size=size(Patt{1}.template,2)-1;
target=movmean(Result_perm.spike,5,2);
for i=1:length(Patt)
    CC=normxcorr2(Patt{i}.template,target);
    xcorr_img(i,:)=CC(size(target,1),:);
end
xcorr_img=xcorr_img(:,wind_size/2+1:end-wind_size/2);
xcorr_img_th=zeros(size(xcorr_img,1),size(xcorr_img,2));
for i=1:size(xcorr_img,1)
    [pks s_tmp width prom]=findpeaks(xcorr_img(i,:),'MinPeakDistance',wind_size/4);
    s_tmp=s_tmp(find(pks>corr_th));
    xcorr_img_th(i,s_tmp)=xcorr_img(i,s_tmp);
end
end