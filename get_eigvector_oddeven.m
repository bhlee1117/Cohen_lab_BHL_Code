function [V, D, eigTrace, validInd, covMat] = get_eigvector_oddeven(subMov1, subMov2, n)
% [V D eigTrace validInd covMat] = get_eigvector_oddeven(subMov1, subMov2, n)
%
% Shot-noise-free PCA: eigendecomposition of the ODD/EVEN CROSS-covariance
% matrix instead of the ordinary covariance matrix.
%
% subMov1 and subMov2 are two measurements of the SAME signal with INDEPENDENT
% noise (e.g. the same ROI/footprint projected onto its odd- and even-checkerboard
% pixel halves, or two interleaved camera channels):
%       x1 = s + n1 ,  x2 = s + n2 ,  <n1 n2'> = 0 ,  <s n_k'> = 0
%
% The ordinary covariance of either channel is
%       C_std = <x x'> = S + N
% where N is diagonal-dominant shot noise; it inflates the diagonal and makes the
% PCA variance spectrum flatten into a noise floor (all late PCs ~ equal variance).
% The cross-covariance instead is
%       C_oe  = <x1 x2'> = S           (in expectation; noise cancels, incl. the diagonal)
% so the eigenvalue spectrum keeps falling instead of hitting a floor, and the
% leading eigenvectors are not contaminated by per-ROI noise.
%
% Because the noise only cancels IN EXPECTATION, finite-T estimates of C_oe are
% not positive semi-definite: small eigenvalues scatter around zero and some come
% out NEGATIVE. That is the expected behaviour -- the magnitude of the negative
% tail is the estimator's noise floor (|D| there ~ sampling error), and any
% eigenvalue clearly above it is real, shot-noise-free signal variance.
%
% INPUT
%   subMov1, subMov2 - [N x T] traces (ROIs/pixels x time), NaN allowed.
%                      Same signal, independent noise.
%   n                - (optional) number of components to return. Default: all.
%
% OUTPUT
%   V        - [N x n]  eigenvectors (spatial modes), sorted by descending D.
%                       Sign convention matches get_eigvector.
%   D        - [N x 1]  eigenvalues (shot-noise-free variance per component).
%                       May be negative at the tail (see above).
%   eigTrace - [Tvalid x n] temporal scores of the COMBINED trace (x1+x2)/2
%                       projected onto V.
%   validInd - [1 x T] logical, frames used (no NaN in either channel).
%   covMat   - [N x N]  the symmetrized odd/even cross-covariance that was
%                       diagonalized.
%
% Validation (synthetic, N=40 ROIs, T=2e4, K=3 shared components):
%   isotropic noise      -> both methods find the same subspace, but the ordinary
%                           eigenvalues are inflated by the noise variance and the
%                           tail plateaus at that floor, hiding where K ends.
%   heteroscedastic noise (bright vs dim footprints, the real case)
%                        -> ordinary PCA subspace error 33 deg, odd/even 3.6 deg;
%                           odd/even also recovers the true eigenvalues to ~1%.
%
% See also: get_eigvector (ordinary covariance version), get_corrTensor
%
% 2026.08.25, Byung Hun Lee

if nargin < 2 || isempty(subMov2)
    error('get_eigvector_oddeven requires two channels with independent noise.');
end
if ~isequal(size(subMov1), size(subMov2))
    error('subMov1 and subMov2 must be the same size ([N x T]).');
end

validInd = sum(isnan(subMov1), 1) == 0 & sum(isnan(subMov2), 1) == 0;
X1 = subMov1(:, validInd);
X2 = subMov2(:, validInd);
T  = size(X1, 2);
if T < 2
    error('Fewer than 2 valid (NaN-free) frames.');
end

X1 = X1 - mean(X1, 2);
X2 = X2 - mean(X2, 2);

% symmetrized odd/even cross-covariance: <x1 x2'> = signal covariance only
covMat = (X1 * X2') / (T - 1);
covMat = (covMat + covMat') / 2;

[V, D] = eig(covMat, 'vector');
[D, si] = sort(D, 'descend');          % eig is not sorted for indefinite matrices
V = V(:, si);

vSign = sign(max(V) - max(-V));        % make the largest value always positive
vSign(vSign == 0) = 1;
V = V .* vSign;

if nargin < 3 || isempty(n)
    n = size(V, 2);
end
n = min(n, size(V, 2));
V = V(:, 1:n);

% temporal scores from the combined (noise-averaged) trace
eigTrace = ((X1 + X2) / 2)' * V;
end
