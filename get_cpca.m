function [U, pc, varFrac, Z] = get_cpca(Y, fs, band, nPC, goodMask)
% [U, pc, varFrac, Z] = get_cpca(Y, fs, band, nPC, goodMask)
%
% COMPLEX PCA: band-limit each pixel's trace, take its analytic signal, then SVD the
% complex matrix.
%
% WHY, FOR A TRAVELLING WAVE. Ordinary PCA on a propagating disturbance cannot put it in
% one component. A wave A*cos(kx - wt) is cos(kx)cos(wt) + sin(kx)sin(wt) -- two separable
% terms -- so real PCA returns TWO modes, spatially 90 deg apart and temporally in
% quadrature, and splits the wave's variance between them. Neither mode alone is the
% wave, the split point depends on noise, and the eigenvalue ordering can swap the two.
%
% The analytic signal removes the degeneracy. Writing each pixel as z(x,t) = a(x)e^{i0(x)}
% * e^{-iwt}, a travelling wave is rank ONE in the complex domain: a single mode whose
% spatial eigenvector carries the amplitude in |U| and the propagation in angle(U). The
% phase gradient along the vessel IS the wave, read directly off the eigenvector, and the
% mode's time course is already the analytic A(t)e^{i phi(t)} -- so real(pc) and imag(pc)
% are a quadrature pair with no second Hilbert transform needed.
%
% BAND-LIMIT FIRST, ALWAYS. The analytic signal of a broadband trace has no meaningful
% instantaneous phase (Bedrosian's condition fails), so `band` is not optional -- pass the
% band the oscillation is expected in.
%
% PHASE REFERENCE. A complex mode is defined only up to a global rotation e^{i theta},
% the complex analogue of PCA's arbitrary sign. Each mode is rotated so its
% amplitude-weighted mean phase is zero, which makes angle(U) comparable between segments
% and across runs. The rotation is applied to U and undone on pc, so U*pc is unchanged.
%
% INPUT
%   Y         [nPix x nT]  real traces, one row per pixel (mean-removed by the caller)
%   fs        sampling rate (Hz)
%   band      [lo hi] Hz, the passband applied before the Hilbert transform
%   nPC       number of modes to return
%   goodMask  logical [1 x nT], frames the SVD is computed on ([] = all). pc is
%             returned over ALL frames regardless, so bad frames are excluded from the
%             mode estimate without punching holes in the time courses.
%
% OUTPUT
%   U         [nPix x nPC]  complex spatial eigenvectors, phase-referenced
%   pc        [nPC x nT]    complex mode time courses over every frame
%   varFrac   [nPC x 1]     fraction of the retained variance in each mode
%   Z         [nPix x nT]   the analytic signal the SVD was taken of
%
% 2026.08.28, Byung Hun Lee

if nargin < 5 || isempty(goodMask), goodMask = true(1, size(Y,2)); end
assert(numel(band) == 2 && band(2) > band(1), 'band must be [lo hi] with hi > lo.');

Yb = get_bandpass(Y, fs, band);
Yb(~isfinite(Yb)) = 0;

Z = hilbert(Yb.').';                       % hilbert works down columns -> transpose
Z = Z - mean(Z, 2);

nPC = max(1, min(nPC, min(size(Z,1), nnz(goodMask)) - 1));
[U, Sv, ~] = svds(Z(:, goodMask), nPC);

pc = U' * Z;                               % ' is the conjugate transpose: correct here

% fix the arbitrary global rotation of each mode
for j = 1:size(U,2)
    th = angle(sum(U(:,j) .* abs(U(:,j))));       % amplitude-weighted mean phase
    if ~isfinite(th), th = 0; end
    U(:,j)  = U(:,j)  * exp(-1i*th);
    pc(j,:) = pc(j,:) * exp( 1i*th);              % keeps U*pc invariant
end

varFrac = diag(Sv).^2;
varFrac = varFrac / max(sum(varFrac), eps);
end
