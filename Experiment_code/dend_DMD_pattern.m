function dend_DMD_pattern(app,thres)

dmd_now=app.GetDevice('DMD_Device');

[fname fpath]=uigetfile();
load(fullfile(fpath,fname))
im=snap.img;
im=medfilt2(im,[5 5],"symmetric");
imfilt=mat2gray(im-medfilt2(im,[150 150],"symmetric"));
figure(2); clf;
nexttile([1 1]); imshow2(im,[])
nexttile([1 1]); imshow2(imfilt,[])

% Apply watershed
%L = watershed(imfilt);
im_filt_th=imfilt-prctile(imfilt(imfilt>0),thres);
im_filt_th(im_filt_th<0)=0;
L = bwlabel(im_filt_th);
[Npixel Nws]=histcounts(L(:),[0:max(L(:))]);
candN=Nws(Npixel>1000);
avgVal=zeros(length(candN),1);
for n=2:length(candN)
    regions=imfilt(L==candN(n));
    avgVal(n)=median(regions);
end
candN=candN(find(avgVal>0.01));
[mask L_mask]=ismember(L,candN);
L_mask(~mask)=0;
L_maskJet = label2rgb(L_mask+1,'turbo',[.5 .5 .5]);
nexttile([1 1]); imshow2(L_maskJet,[])

L_mask_filt=zeros(size(L_mask));
for m=1:max(L_mask(:))
    L_mask_filt(imgaussfilt(double(L_mask==m),8) > 0.2)=1;
end
nexttile([1 1]); imshow2(L_mask_filt,[])
end


