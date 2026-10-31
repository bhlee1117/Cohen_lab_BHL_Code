function [C, R, nPair] = get_lag1cov(X, mask)
% get_lag1cov  Lag-1 covariance and correlation between rows of X.
%
%   [C, R, nPair] = get_lag1cov(X, mask)
%
% Estimates the ROI x ROI covariance from the ONE-FRAME LAG product
% X(t)X(t+1)' instead of the same-time product X(t)X(t)'. Photon shot noise is
% temporally white, so it contributes nothing to the lag-1 product -- including
% on the DIAGONAL, which is where same-time covariance puts all of it. That is
% what makes this an alternative to the odd/even (cross-channel) trick: it needs
% only one trace per ROI, so it works on traces that cannot be split into two
% independent-noise halves.
%
% IMPORTANT: this only cancels shot noise on UNSMOOTHED traces. After a k-frame
% moving average the noise has lag-1 autocorrelation ~ (k-1)/k, so most of it
% survives into the lag-1 product. Pass spike-removed but unsmoothed traces
% (e.g. get_subthreshold(..., 1)) rather than 20-ms subthreshold traces.
%
% The price is that real signal decorrelating within one frame is attenuated by
% roughly exp(-1/tau). Slow subthreshold signal is essentially unaffected.
%
% Inputs
%   X    : N x T traces (one row per ROI), NaN allowed.
%   mask : 1 x T logical, frames to analyse. Default: all frames.
%          Frames outside the mask break lag-1 pairs rather than being closed
%          up, so a lag of one frame is always a true one-frame separation.
%
% Outputs
%   C     : N x N lag-1 covariance, symmetrised as (A1*A2' + A2*A1')/2.
%   R     : N x N correlation, C normalised by its own diagonal. Rows/columns
%           whose lag-1 variance is <= 0 come back NaN -- that is the unbiased
%           estimator landing below zero for an ROI with no reliable signal,
%           not a failure.
%   nPair : number of consecutive analysable frame pairs used.
%
% See also: get_corrMat, get_corrTensor, covEstimators3

if nargin < 2 || isempty(mask), mask = true(1, size(X,2)); end
mask = logical(mask(:)');
assert(numel(mask) == size(X,2), 'get_lag1cov:mask', 'mask must be 1 x T');

g = mask & all(isfinite(X), 1);            % one common frame set for every ROI
pr = find(g(1:end-1) & g(2:end));          % both frames of the pair analysable
nPair = numel(pr);
N = size(X,1);
assert(nPair > N + 2, 'get_lag1cov:few', ...
       'only %d consecutive analysable frame pairs for %d ROIs', nPair, N);

Xm = X - mean(X(:, g), 2);                 % ONE mean, over all analysed frames
A1 = Xm(:, pr);  A2 = Xm(:, pr+1);
C  = (A1*A2' + A2*A1') / 2 / (nPair - 1);  % X(t)X(t+1)' is not symmetric
C  = (C + C')/2;                           % kill the last bit of round-off

d  = diag(C);
sd = sqrt(d);  sd(d <= 0) = NaN;           % unusable rows -> NaN, not complex
R  = C ./ (sd * sd');
end
