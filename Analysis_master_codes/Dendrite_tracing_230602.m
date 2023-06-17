load('/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/20230519_DoubleFocus_VU/Snaps/105530test2_somaCA1.mat')
im=snap.img;
im=medfilt2(im,[5 5],"symmetric");
imfilt=mat2gray(im-medfilt2(im,[150 150],"symmetric"));
im=mat2gray(im);
figure(2); clf;
nexttile([1 1]); imshow2(im,[])
nexttile([1 1]); imshow2(imfilt,[])

% Apply watershed
%L = watershed(imfilt);
im_filt_th=imfilt-prctile(imfilt(imfilt>0),50);
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



%%
FOV_MM = 1; % assume square FOV
% Parameters to set

OPT.MIN_DIAM_UM = 3;
OPT.MAX_DIAM_UM = 8;
OPT.FRANGI_SCALE_RATIO = 1;
OPT.FRANGI_BETA_ONE = 0.5;
OPT.FRANGI_BETA_TWO = 12;
OPT.BINARIZATION_LEVEL = 0.2;

MM_per_PU = FOV_MM ./ size(imfilt);
min_vessel_diameter_PU = (OPT.MIN_DIAM_UM / MM_per_PU(1)) / 1E3;
max_vessel_diameter_PU = (OPT.MAX_DIAM_UM / MM_per_PU(1)) / 1E3;
min_s = ceil(min_vessel_diameter_PU / 2);
max_s = ceil(max_vessel_diameter_PU / 2);

options.BlackWhite = false;
options.FrangiScaleRatio = OPT.FRANGI_SCALE_RATIO;
options.FrangiScaleRange = [min_s max_s];
options.FrangiBetaOne = OPT.FRANGI_BETA_ONE;
options.FrangiBetaTwo = OPT.FRANGI_BETA_TWO;

[outIm,~,~] = FrangiFilter2D(255*imfilt, options);
large_vessel_mask = outIm > OPT.BINARIZATION_LEVEL;
se = strel('disk', 5); % Define a structuring element
skeletonImage = bwmorph(large_vessel_mask, 'skel', Inf);

figure; clf;
ax1=[];
ax1=[ax1 nexttile([1 1])];
imshow2(labeloverlay(imfilt,skeletonImage,'Transparency',0),[])
ax1=[ax1 nexttile([1 1])];
imshow2(large_vessel_mask,[])
ax1=[ax1 nexttile([1 1])];
imshow2(imfuse(large_vessel_mask,imfilt),[])
ax1=[ax1 nexttile([1 1])];
imshow2(outIm,[])
linkaxes(ax1,'xy');

%%
figure(1); clf;
MergeIm=outIm.*imfilt;
%MergeIm(MergeIm<0.001)=0;
[~,threshold] = edge(MergeIm,'sobel');
fudgeFactor = 0.5;
BWs = edge(MergeIm,'sobel',threshold * fudgeFactor);
ax1=[];
ax1=[ax1 nexttile([1 1])];
imshow2(labeloverlay(outIm,BWs,'Transparency',0),[])
se90 = strel('diamond',2);
se0 = strel('line',3,0);

BWsdil = imdilate(BWs,[se90 se0]);
b=bwlabel(BWsdil);
for i = 1:max(b(:))
ind=find(b==i);
Area(i)=length(ind);
point(i)=ind(1);
end
Bwinterest=find(Area>15000);
[pointY pointX]=ind2sub(size(outIm),point(Bwinterest)); 
ax1=[ax1 nexttile([1 1])];
imshow2(labeloverlay(im,BWsdil,'Transparency',0),[])
hold all
plot(pointX, pointY,'ro')

%BWdfill = imfill(BWsdil,point(Bwinterest)');
%BWdfill = imfill(BWsdil);
BWdfill=BWsdil;
BWnobord = imclearborder(BWdfill,1);

ax1=[ax1 nexttile([1 1])];
imshow2(labeloverlay(im,BWdfill,'Transparency',0),[])



seD = strel('diamond',4);
BWfinal = imerode(BWnobord,seD);
BWfinal = imerode(BWfinal,seD);
BWfinal = imdilate(BWfinal,[se90 se0]);
ax1=[ax1 nexttile([1 1])];
imshow2(labeloverlay(outIm,BWfinal,'Transparency',0),[])
linkaxes(ax1,'xy');