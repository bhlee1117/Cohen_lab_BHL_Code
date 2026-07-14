function [crossCorr, lags, nPairs] = nanXCorrFFT(x, y, maxLag, normalized, minPairs)
% nanXCorrFFT  FFT-based NaN-aware cross-correlation (fast nanXCorr).
%
% Identical result and convention to nanXCorr, but O(n log n) instead of
% O(n * maxLag): every lag is obtained from a handful of FFT cross-
% correlations rather than a per-lag loop.
%
% Convention (same as nanXCorr):
%   positive lag k means y is delayed relative to x:
%       crossCorr(k) = corr( x(t), y(t+k) )   (normalized = true)
%       crossCorr(k) = mean( x(t) .* y(t+k) ) (normalized = false)
%   with means/variances computed per lag over the valid (non-NaN) overlap.
%
% INPUT / OUTPUT identical to nanXCorr:
%   [crossCorr, lags, nPairs] = nanXCorrFFT(x, y, maxLag, normalized, minPairs)
%
% See also: nanXCorr

x = x(:);  y = y(:);
n = numel(x);
if numel(y) ~= n, error('Input vectors must be of the same length.'); end
if n == 0, error('Input vectors must not be empty.'); end
if nargin < 3 || isempty(maxLag),     maxLag = n - 1;  end
if nargin < 4 || isempty(normalized), normalized = true; end
if nargin < 5 || isempty(minPairs),   minPairs = 3;    end
maxLag = min(maxLag, n - 1);

% Validity masks and NaN-zeroed signals
mx = double(~isnan(x));   xz = x;  xz(mx == 0) = 0;
my = double(~isnan(y));   yz = y;  yz(my == 0) = 0;

% Zero-pad so circular correlation equals linear correlation for |k| <= maxLag
L = 2^nextpow2(n + maxLag);

% FFTs: conjugate on the x-side (left) operands, plain on the y-side (right)
Fxz  = conj(fft(xz,     L));   Gyz  = fft(yz,     L);
Fxz2 = conj(fft(xz.^2,  L));   Gyz2 = fft(yz.^2,  L);
Fmx  = conj(fft(mx,     L));   Gmy  = fft(my,     L);

% Reorder a length-L circular correlation into lags -maxLag..maxLag.
% index m (0-based) of ifft(conj(FA).*GB) = sum_t A(t) B(t+m).
pick = @(c) [c(L-maxLag+1:L); c(1:maxLag+1)];
xc   = @(FA, GB) pick(real(ifft(FA .* GB)));

% Per-lag sums over the valid overlap
nP  = round(xc(Fmx,  Gmy));    % # valid pairs
Sxy = xc(Fxz,  Gyz);           % sum x(t) y(t+k)
Sx  = xc(Fxz,  Gmy);           % sum x(t) over valid pairs
Sy  = xc(Fmx,  Gyz);           % sum y(t+k) over valid pairs
Sxx = xc(Fxz2, Gmy);           % sum x(t)^2 over valid pairs
Syy = xc(Fmx,  Gyz2);          % sum y(t+k)^2 over valid pairs

lags   = (-maxLag:maxLag).';
nPairs = nP;

if normalized
    covXY = Sxy - (Sx .* Sy) ./ nP;
    varX  = Sxx - (Sx.^2)    ./ nP;
    varY  = Syy - (Sy.^2)    ./ nP;
    crossCorr = covXY ./ sqrt(varX .* varY);
    crossCorr(varX <= 0 | varY <= 0) = NaN;
else
    crossCorr = Sxy ./ nP;     % mean product
end

crossCorr(nP < minPairs) = NaN;
end
