function [Cmat, lags] = get_corrTensor(tr1, tr2, trC, mask, maxLag, minPairs)
% get_corrTensor  Pairwise spatiotemporal correlation tensor (N x N x (2*maxLag+1)).
%
%   [Cmat, lags] = get_corrTensor(tr1, tr2, trC, mask, maxLag, minPairs)
%
% Cmat(i,j,:) is the cross-correlogram between ROI i and ROI j vs lag, computed
% on the RAW traces over the frames selected by MASK. Non-selected frames are
% set to NaN so the time axis stays continuous (a lag of k is a true k-frame
% separation, no concatenation across gaps). Same-ROI pairs (diagonal) cross the
% two optical channels tr1 x tr2 to cancel the zero-lag shot-noise spike;
% different ROIs use the combined trace trC (their noise is already independent).
% Correlations use nanXCorrFFT.
%
%   * spatial correlation matrix = zero-lag slice  Cmat(:,:,maxLag+1)
%   * temporal correlation       = the lag profile Cmat(i,j,:)
%
% Inputs:
%   tr1, tr2  - N x T per-channel traces (diagonal / autocorrelation)
%   trC       - N x T combined trace     (off-diagonal / different ROIs)
%   mask      - 1 x T logical, frames to include
%   maxLag    - max lag in frames
%   minPairs  - min valid frame-pairs per lag (passed to nanXCorrFFT)
%
% See also: nanXCorrFFT, get_corrMatAutoCh

nROI = size(trC, 1);
nLag = 2*maxLag + 1;
mask = logical(mask(:)');
tr1(:, ~mask) = NaN;  tr2(:, ~mask) = NaN;  trC(:, ~mask) = NaN;

Cmat = nan(nROI, nROI, nLag);
lags = (-maxLag:maxLag).';
for i = 1:nROI
    for j = i:nROI
        if i == j              % same ROI -> two channels (noise-free at zero lag)
            a = tr1(i, :);  b = tr2(j, :);
        else                   % different ROIs -> combined trace
            a = trC(i, :);  b = trC(j, :);
        end
        if sum(~isnan(a) & ~isnan(b)) <= 2*maxLag, continue; end
        xc = nanXCorrFFT(a, b, maxLag, true, minPairs);
        Cmat(i, j, :) = xc;
        if i ~= j, Cmat(j, i, :) = flip(xc); end   % C(j,i,tau) = C(i,j,-tau)
    end
end
end
