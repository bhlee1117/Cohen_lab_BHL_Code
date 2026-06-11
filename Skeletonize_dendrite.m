function DMDmask = Skeletonize_dendrite(im, thickness, thres, mag, showFig)
% SKELETONIZE_DENDRITE  Detect and skeletonize dendrite/vessel structures
%                       in a fluorescence image using Frangi vesselness filtering.
%
% Usage:
%   DMDmask = Skeletonize_dendrite(im, thickness, thres)
%   DMDmask = Skeletonize_dendrite(im, thickness, thres, mag)
%   DMDmask = Skeletonize_dendrite(im, thickness, thres, mag, showFig)
%
% Description:
%   Applies a difference-of-Gaussians (DoG) pre-filter to enhance tubular
%   structures, then runs a 2D Frangi vesselness filter to detect dendrite-
%   like structures. The detected mask is skeletonized and dilated to the
%   specified thickness, producing a binary mask suitable for DMD targeting.
%
% Inputs:
%   im        - 2D grayscale input image
%   thickness - Dilation thickness (in pixels) applied to the skeleton
%   thres     - Binarization threshold for the Frangi filter output (0–1)
%   mag       - Objective magnification (default: 25)
%   showFig   - Logical flag to display the result figure (default: true)
%
% Output:
%   DMDmask   - Binary mask of the dilated dendrite skeleton

if nargin < 4
    mag = 25;
end
if nargin < 5
    showFig = true;
end

% --- Edge boundary mask (suppress edge artifacts) ---
bound = 6;
bound_mask = ones(size(im, 1), size(im, 2));
bound_mask(bound:end-bound, bound:end-bound) = 0;

% --- Pre-filtering: Gaussian smoothing + DoG enhancement ---
im = imgaussfilt(im, 0.2 * mag / 10);
imfilt = mat2gray(imgaussfilt(im - imgaussfilt(im, 20), 0.7));
im = mat2gray(im);

% --- Frangi vesselness filter parameters ---
FOV_MM = 0.5;                       % Field of view in mm (assumes 0.5 mm FOV)
OPT.MIN_DIAM_UM        = 0.1;       % Minimum vessel diameter (um)
OPT.MAX_DIAM_UM        = 6;         % Maximum vessel diameter (um)
OPT.FRANGI_SCALE_RATIO = 1;
OPT.FRANGI_BETA_ONE    = 0.5;       % Sensitivity to blob vs. line structures
OPT.FRANGI_BETA_TWO    = 12;        % Sensitivity to background noise
OPT.BINARIZATION_LEVEL = thres;

MM_per_PU = FOV_MM ./ size(imfilt);
min_vessel_diameter_PU = (OPT.MIN_DIAM_UM / MM_per_PU(1)) / 1E3;
max_vessel_diameter_PU = (OPT.MAX_DIAM_UM / MM_per_PU(1)) / 1E3;
min_s = ceil(min_vessel_diameter_PU / 2);
max_s = ceil(max_vessel_diameter_PU / 2);

options.BlackWhite        = false;
options.FrangiScaleRatio  = OPT.FRANGI_SCALE_RATIO;
options.FrangiScaleRange  = [min_s max_s];
options.FrangiBetaOne     = OPT.FRANGI_BETA_ONE;
options.FrangiBetaTwo     = OPT.FRANGI_BETA_TWO;

% --- Frangi filtering and binarization ---
[outIm, ~, ~] = FrangiFilter2D(255 * imfilt, options);
large_vessel_mask = outIm > OPT.BINARIZATION_LEVEL;
large_vessel_mask(bound_mask == 1) = 0;

% --- Skeletonize and dilate ---
skeletonImage = bwmorph(large_vessel_mask, 'skel', Inf);

se90   = strel('line', thickness, 90);
se0    = strel('line', thickness, 0);
DMDmask = imdilate(skeletonImage, [se90 se0]);

% --- Optional display ---
if showFig
    figure;
    imshow2(imfuse(im, DMDmask), []);
    title('Skeletonized Dendrite Mask');
end

end
