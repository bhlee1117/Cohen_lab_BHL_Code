function [dtimg] = get_dtimg(STAmovie)
%get_dtimg  Generate spatiotemporal SNAPt movie with subframe timing.
%
%   [dtimg] = get_dtimg(STAmovie)
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
%
%   OUTPUTS:
%     dtimg_Reg        – Registered spike timing map (Y × X).
%     dFF_Reg          – Registered ΔF/F image (Y × X).
%
%   Byung-Hun, 2025.

% === Pre-processing ===
spikeMov  = STAmovie;
dSpikeMov = spikeMov - mean(spikeMov,3);              % deviations from mean

% === Spatial & PCA filtering ===
filtsize = 7; sigmaSpatial = 2;
spikeFilt = spatialfilt(dSpikeMov, filtsize, sigmaSpatial);
spikeFilt(isnan(spikeFilt)) = 0;
[spikeFiltPCA, eigvecs, eigvals] = pcafilt(spikeFilt, [1:5]);

% === Spline fit to estimate timing ===
disp('Fitting spline to estimate spike timing (computationally demanding)…');
tic
thresh = 0.5; dir = 1;
dtimg = splinezeros(spikeFiltPCA, thresh, dir);
dtimg = dtimg - nanmean(dtimg(:));
dtimg = set_edge(dtimg,2,NaN);
validImg = imfill(~isnan(dtimg),'holes');
toc

% === Display timing & amplitude maps ===
figure; clf;
nexttile([1 1])
dtimg(dtimg==0)=NaN;
imshow2(dtimg - prctile(dtimg(:),5), [0 3]); title('Timing')
end
