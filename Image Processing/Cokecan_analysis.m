function [intensity1] = Cokecan_analysis(mov1);

Initialblur = 0.5;  % blur to apply before doing erosion to define background
bkgdisksize = 7;    % radius of the disk used for erosion to define background.  Should be just bigger than one cell.
Watershedblur = 1;  % blur to apply before watershed
Regionthresh = 0.5;  % fraction of max - min intensity for defining a mask within a region


% Open the experimental parameters file and get the image size
%Info=textscan(fopen('experimental_parameters.txt'),'%s');


firstpic=mov1(:,:,1);
npix = 5;       % smoothing function is an npix by npix square Gaussian weighting with a stdev of sigma
sigma = Watershedblur;
GaussFilt = fspecial('gaussian', npix, sigma);
pic3 = imfilter(firstpic,GaussFilt,'replicate');    
imshow(pic3,[])


L = watershed(-pic3, 8);  % Find the watershed regions in the image
figure(5)
pcolor(double(L))
colormap('jet')
shading 'interp'
title('Watershed regions')

stats = regionprops(L, firstpic, 'Area', 'MaxIntensity', 'MinIntensity', 'PixelIdxList');

hist([stats(:).Area], 100)  % histogram the areas of the watershed regions
set(gca, 'yscale', 'log')

toosmallthresh = input('Input the minimum area to keep: ');

toosmallindx = find([stats(:).Area] < toosmallthresh); 

stats2 = stats;
stats2(toosmallindx) = [];
figure(6)
hist([stats2(:).MaxIntensity], 100);
set(gca, 'yscale', 'log')

toodimthresh = input('Input the minimum intensity to keep: ');

toodimindx = find([stats2(:).MaxIntensity] < toodimthresh);
stats2(toodimindx) = [];

figure(7)
plot([stats2(:).MaxIntensity], [stats2(:).Area], '.')
xlabel('Maximum intensity in each region')
ylabel('Region area')

% Threshold each region based on its local contrast.
ncells = length(stats2);
L2 = zeros(size(L));
for j = 1:ncells;
    maxin = stats2(j).MaxIntensity;
    minin = stats2(j).MinIntensity;
    thresh = Regionthresh*(maxin - minin) + minin;
    abovethreshIdx = pic3(stats2(j).PixelIdxList) > thresh;
    L2(stats2(j).PixelIdxList(abovethreshIdx)) = j;
end;
imshow(L2);
binL2 = sign(L2);
se = strel('disk', 1);
binL3 = imclose(binL2, se); % join regions that are too close together

imshow(binL3)
L4 = bwlabel(binL3);
pcolor(L4); shading 'interp'; colormap 'jet'

stats4 = regionprops(L4, firstpic, 'Area', 'MeanIntensity', 'PixelIdxList');


% Select the region of cells to keep
figure(4); clf
subplot(1,2,1);
plot([stats4(:).Area], [stats4(:).MeanIntensity], '.');
xlabel('Area'); ylabel('Mean Intensity');
title('Select the region to keep, right click to end')
subplot(1,2,2);
imshow(firstpic, []);
    subplot(1,2,1)
    [xvnew, yvnew] = getline(gcf, 'closed');
  if length(xvnew) < 3
        return;
    else
        xv = xvnew;
        yv = yvnew;
    end;
    inpoint = inpolygon([stats4(:).Area], [stats4(:).MeanIntensity], xv, yv);

    stats5 = stats4(inpoint);

    L5 = zeros(size(L4));
    nkeep = length(stats5);
    for j = 1:nkeep
        L5(stats5(j).PixelIdxList) = 1;
    end;

    colorimg = repmat(mat2gray(firstpic), [1 1 3]);
    colorimg(:,:,2) = colorimg(:,:,2) + .5*L5;
    colorimg(:,:,3) = colorimg(:,:,3).*(1- L5);
    colorimg(:,:,1) = colorimg(:,:,1).*(1- L5);
    subplot(1,2,2)
    imshow(colorimg, []);

imwrite(sign(L5), ['mask_' datestr(now, 30) '.tif'], 'tif');

cc=bwconncomp(L5); % select ROI 
intensity1=zeros(cc.NumObjects,size(mov1,3));

for k=1:size(mov1,3)
x=regionprops(cc,mov1(:,:,k),'MeanIntensity');
intensity1(:,k)=cell2mat(struct2cell(x));
end





