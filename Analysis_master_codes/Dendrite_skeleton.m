FOV_MM = 1; % assume square FOV
% Parameters to set
im=double(im);
im=medfilt2(im,[5 5],"symmetric");
imfilt=mat2gray(im-medfilt2(im,[150 150],"symmetric"));

OPT.VESSEL_REMOVAL_MIN_DIAM_UM = 5;
OPT.VESSEL_REMOVAL_MAX_DIAM_UM = 10;
OPT.VESSEL_REMOVAL_FRANGI_SCALE_RATIO = 1;
OPT.VESSEL_REMOVAL_FRANGI_BETA_ONE = 0.5;
OPT.VESSEL_REMOVAL_FRANGI_BETA_TWO = 12;
OPT.VESSEL_REMOVAL_BINARIZATION_LEVEL = 0.2;

MM_per_PU = FOV_MM ./ size(imfilt);
min_vessel_diameter_PU = (OPT.VESSEL_REMOVAL_MIN_DIAM_UM / MM_per_PU(1)) / 1E3;
max_vessel_diameter_PU = (OPT.VESSEL_REMOVAL_MAX_DIAM_UM / MM_per_PU(1)) / 1E3;
min_s = ceil(min_vessel_diameter_PU / 2);
max_s = ceil(max_vessel_diameter_PU / 2);

options.BlackWhite = false;
options.FrangiScaleRatio = OPT.VESSEL_REMOVAL_FRANGI_SCALE_RATIO;
options.FrangiScaleRange = [min_s max_s];
options.FrangiBetaOne = OPT.VESSEL_REMOVAL_FRANGI_BETA_ONE;
options.FrangiBetaTwo = OPT.VESSEL_REMOVAL_FRANGI_BETA_TWO;

[outIm,~,~] = FrangiFilter2D(255*imfilt, options);
large_vessel_mask = outIm > OPT.VESSEL_REMOVAL_BINARIZATION_LEVEL;
figure(1); clf;
ax1=[];
ax1=[ax1 nexttile([1 1])];
imshow2(imfilt,[])
ax1=[ax1 nexttile([1 1])];
imshow2(large_vessel_mask,[])
ax1=[ax1 nexttile([1 1])];
imshow2(imfuse(large_vessel_mask,imfilt),[])
ax1=[ax1 nexttile([1 1])];
imshow2(outIm,[])
linkaxes(ax1,'xy');

%%
figure;
se = strel('disk', 5); % Define a structuring element
%cleanImage = imopen(outIm, se); % Perform morphological opening
skeletonImage = bwmorph(large_vessel_mask, 'skel', Inf);
imshow2(skeletonImage, []);
hold off;