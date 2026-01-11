function [superlocColormov, dtimg_Reg, dFF_Reg] = generate_SNAPTmov(STAmovie, mask, StrImg, tformReg)
%GENERATE_SNAPTMOV  Generate spatiotemporal SNAPt movie with subframe timing.
%
%   [superlocColormov, dtimg_Reg, dFF_Reg] = GENERATE_SNAPTmov(STAmovie, mask, StrImg, tformReg)
%
%   This function creates a super-resolved temporal (SNAPT) movie of
%   action-potential timing across a masked ROI. It:
%     • Computes ΔF/F image and spike-triggered deviations
%     • Applies spatial + PCA filtering to suppress noise
%     • Fits a spline to estimate per-pixel spike timing (dtimg)
%     • Optionally applies an affine/tform registration to align images
%     • Generates a time-resolved Gaussian "flash" movie weighted by
%       standard intensity and ΔF/F
%     • Overlays the coloured SNAPT movie on the structural image
%
%   INPUTS:
%     STAmovie  – 3-D array (Y × X × T), spike-triggered average movie.
%     mask      – Binary mask of ROI (Y × X). Single plane only.
%     StrImg    – Structural/brightfield image for background (Y × X).
%     tformReg  – (Optional) geometric transform (e.g. from imregtform)
%
%   OUTPUTS:
%     superlocColormov – 4-D RGB movie (Y × X × 3 × Nframes).
%     dtimg_Reg        – Registered spike timing map (Y × X).
%     dFF_Reg          – Registered ΔF/F image (Y × X).
%
%   Byung-Hun's lab, 2025.

if nargin < 4
    tformReg = [];
end
if size(mask,3) > 1
    error('Mask must be a single plane (Y×X).');
end

% === Pre-processing ===
spikeMov  = STAmovie;
dSpikeMov = spikeMov - mean(spikeMov,3);              % deviations from mean
bkgImg    = mean(STAmovie(:,:,end-5:end),3);          % background image
dFMov     = STAmovie - bkgImg;
APamp     = imgaussfilt(max(dFMov,[],3),5);
dFF       = max(APamp ./ imgaussfilt(bkgImg,5), [],3);

% === Spatial & PCA filtering ===
filtsize = 7; sigmaSpatial = 2;
spikeFilt = spatialfilt(dSpikeMov, filtsize, sigmaSpatial);
spikeFilt(isnan(spikeFilt)) = 0;
[spikeFiltPCA, eigvecs, eigvals] = pcafilt(spikeFilt, 5);

% === Spline fit to estimate timing ===
disp('Fitting spline to estimate spike timing (computationally demanding)…');
tic
thresh = 0.1; dir = 1;
dtimg = splinezeros(spikeFiltPCA, thresh, dir);
dtimg = dtimg - nanmean(dtimg(:));
dtimg = set_edge(dtimg,5,NaN);
validImg = imfill(~isnan(dtimg),'holes');
toc

% === Registration (if provided) ===
if ~isempty(tformReg)
    ref = imref2d(size(StrImg));
    dtimg_Reg   = imwarp(dtimg,       tformReg, 'OutputView', ref);
    spikeMovReg = imwarp(dSpikeMov,   tformReg, 'OutputView', ref);
    dFF_Reg     = imwarp(dFF,         tformReg, 'OutputView', ref);
    validImgReg = double(imwarp(validImg, tformReg, 'OutputView', ref));
    validImgReg(validImgReg==0)=NaN;
    dtimg_Reg   = dtimg_Reg .* validImgReg;

    dFF_Reg = dFF_Reg .* mask;
    dFF_Reg(dFF_Reg==0)=NaN;
else
    dtimg_Reg   = dtimg;
    spikeMovReg = dSpikeMov;
    dFF_Reg     = dFF;
end

% === Display timing & amplitude maps ===
figure; clf;
nexttile([1 1])
dtimg_Reg = dtimg_Reg .* mask;
dtimg_Reg(dtimg_Reg==0)=NaN;
imshow2(dtimg_Reg - prctile(dtimg_Reg(:),5), [0 3]); title('Timing')
hold on
Maskboundary = cell2mat(bwboundaries(mask));
plot(Maskboundary(:,2),Maskboundary(:,1),'r.')
colormap('turbo')

nexttile([1 1])
imshow2(dFF_Reg,[]); hold on
plot(Maskboundary(:,2),Maskboundary(:,1),'r.')
colormap('turbo')

% === Build SNAPT Gaussian-flash movie ===
subframeT = 0.02; % ms per subframe
initialT  = -2;   % ms
finalT    =  2;   % ms
sigmaT    = 0.05; % ms flash width

[ysize, xsize, ~] = size(spikeMovReg);
times   = initialT:subframeT:finalT;
nFrames = numel(times);

stdimg = mean(spikeMovReg,3);
stdimgNorm = mat2gray(stdimg);
dFFNorm    = mat2gray(dFF_Reg);

GaussPeaksmov = zeros(ysize,xsize,nFrames);
for q = 1:nFrames
    GaussPeaksmov(:,:,q) = exp(-(dtimg_Reg-times(q)).^2/(2*sigmaT^2)) .* double(mask);
end

superlocmov = GaussPeaksmov .* repmat(stdimgNorm,[1 1 nFrames]) .* dFFNorm;
superlocmov(isnan(superlocmov)) = 0;

% === Convert to RGB overlay on structural image ===
superlocColormov = zeros(ysize,xsize,3,nFrames);
StrImgGray = grs2rgb(double(StrImg),colormap(gray));
for j = 1:nFrames
    ColorLayer = grs2rgb(double(superlocmov(:,:,j)*4), colormap("hot"));
    superlocColormov(:,:,:,j) = StrImgGray + ColorLayer;
end
end
