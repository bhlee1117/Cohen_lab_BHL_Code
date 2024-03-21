function DMDmask = mask_dendrite(im,thickness,method,thres)

im=imgaussfilt(im,0.3);
sz=size(im);
imfilt=mat2gray(im-medfilt2(im,[40 40],"symmetric"));
%imfilt=mat2gray(imgaussfilt(im-imgaussfilt(im,30),0.7));
im=mat2gray(im);

FOV_MM=1;
OPT.MIN_DIAM_UM = 1;
OPT.MAX_DIAM_UM = 6;
OPT.FRANGI_SCALE_RATIO = 1;
OPT.FRANGI_BETA_ONE = 0.5;
OPT.FRANGI_BETA_TWO = 12;
OPT.BINARIZATION_LEVEL = thres;

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

switch method
    case 'Skeleton'

skeletonImage = bwmorph(large_vessel_mask, 'skel', Inf);
se90 = strel('line',thickness,90);
se0 = strel('line',thickness,0);
BWfinal = imdilate(skeletonImage,[se90 se0]);

    case 'Edge'

figure(1); clf;
MergeIm=outIm.*imfilt;
[~,threshold] = edge(MergeIm,'sobel');
fudgeFactor = 0.5;
BWs = edge(MergeIm,'sobel',threshold * fudgeFactor);
se90 = strel('line',6,90);
se0 = strel('line',6,0);

BWsdil = imdilate(BWs,[se90 se0]);
BWnobord = imclearborder(BWsdil,1);

seD = strel('diamond',4);
BWfinal = imerode(BWnobord,seD);
BWfinal = imerode(BWfinal,seD);

se90 = strel('line',thickness,90);
se0 = strel('line',thickness,0);
BWfinal = imdilate(BWfinal,[se90 se0]);

end
figure(1); clf;
ax1=nexttile([1 1]);
imshow2(outIm,[])
ax2=nexttile([1 1]);
imshow2(labeloverlay(im,BWfinal,'Transparency',0),[])
linkaxes([ax1 ax2],'xy')

roi=drawpolygon(gca);
roi=round(roi.Position); roi(end+1,:)=roi(1,:);
ROImage = poly2mask(roi(:,1)', roi(:,2)', sz(1), sz(2));

figure(1); clf;
ax1=nexttile([1 1]);
imshow2(outIm,[])
ax2=nexttile([1 1]);
imshow2(labeloverlay(im,BWfinal,'Transparency',0),[])
linkaxes([ax1 ax2],'xy')

roi=drawpolygon(gca);
roi=round(roi.Position); roi(end+1,:)=roi(1,:);
ROImage_rmv = poly2mask(roi(:,1)', roi(:,2)', sz(1), sz(2));

DMDmask= (BWfinal | ROImage) & ~ROImage_rmv;

end