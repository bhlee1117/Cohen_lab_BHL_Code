function DMDmask = Skeletonize_dendrite(im,thickness,thres,mag)
if nargin<5
    mag=25;
end
bound=6;
bound_mask=ones(size(im,1),size(im,2));
bound_mask(bound:end-bound,bound:end-bound)=0;

im=imgaussfilt(im,0.2*mag/10);
sz=size(im);
%imfilt=mat2gray(im-medfilt2(im,[30 30]*mag/10,"symmetric"));
imfilt=mat2gray(imgaussfilt(im-imgaussfilt(im,20),0.7));
%imshow2(imfilt,[])
im=mat2gray(im);

FOV_MM=0.5;
OPT.MIN_DIAM_UM = 0.1;
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
large_vessel_mask(bound_mask==1)=0;
skeletonImage = bwmorph(large_vessel_mask, 'skel', Inf);

se90 = strel('line',thickness,90);
se0 = strel('line',thickness,0);
BWfinal = imdilate(skeletonImage,[se90 se0]);

DMDmask= BWfinal;
figure; imshow2(imfuse(im,DMDmask),[]);
end