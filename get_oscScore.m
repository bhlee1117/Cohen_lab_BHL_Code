function [score, f0, detail] = get_oscScore(X, fs, band)
% [score, f0, detail] = get_oscScore(X, fs, band)
%
% How much does each row of X look like an OSCILLATION (as opposed to noise or drift)?
%
% Two independent properties are required of a real oscillation, and the score is their
% product, so a trace has to satisfy BOTH:
%
%   conc  spectral concentration.  Fraction of the IN-BAND power that sits within +/-20%
%         of the spectral peak. A sinusoid -> 1, a flat spectrum -> small.
%   acPk  autocorrelation side lobe.  Height of the largest autocorrelation peak beyond
%         the first zero crossing, computed on the BAND-PASSED trace. A periodic signal
%         repeats, so it has a strong secondary peak; noise does not.
%
%   score = conc * acPk        in [0, 1]
%
% WHY BOTH, AND WHY BAND-LIMITED: each measure alone has a blind spot, and both are
% fooled by slow drift. Measured on synthetic traces (fs = 1 kHz, band 5-200 Hz):
%
%       signal              conc    acPk   score
%       pure 40 Hz sine     1.000   0.994  0.994
%       narrowband 35-45    1.000   0.927  0.927
%       1/f + 40 Hz         0.744   0.783  0.583
%       40 Hz + 1x noise    0.604   0.588  0.355
%       FM 40 Hz +/-40%     0.505   0.410  0.207
%       1/f (drift only)    0.359   0.166  0.059   <- acPk alone scores this 0.999
%       40 Hz + 3x noise    0.218   0.181  0.039
%       white noise         0.286   0.095  0.027   <- conc alone scores this 0.215
%       broadband 10-200    0.038   0.054  0.002
%
% Oscillatory traces land above ~0.2 and non-oscillatory ones below ~0.06, so a
% threshold of 0.1-0.15 separates them. Restricting the peak search to `band` and
% normalising conc by in-band power is what stops pure drift from scoring highly.
%
% INPUT
%   X     [N x T] traces (rows). NaN allowed: NaN frames are dropped per row.
%   fs    sampling rate (Hz)
%   band  [lo hi] Hz, the frequency range an oscillation is expected in.
%         Default [5, fs/2-1].
%
% OUTPUT
%   score  [N x 1] oscillation score in [0,1]
%   f0     [N x 1] spectral peak inside `band` (Hz)
%   detail struct with fields conc and acPk ([N x 1] each)
%
% See also: get_bandpass, pwelch, xcorr
%
% 2026.08.27, Byung Hun Lee

if nargin < 3 || isempty(band), band = [5, fs/2 - 1]; end
band = [max(band(1), eps), min(band(2), fs/2 - 1)];
assert(band(2) > band(1), 'band must be increasing and below Nyquist.');

if isvector(X), X = X(:).'; end
[N, T] = size(X);
score = nan(N,1);  f0 = nan(N,1);  conc = nan(N,1);  acPk = nan(N,1);

nWin  = min(1024, 2^floor(log2(max(T/2, 8))));
maxLg = min(round(3*fs/band(1)), floor(T/3));      % ~3 cycles of the slowest frequency

for n = 1:N
    x = X(n, :);
    x = x(isfinite(x));
    if numel(x) < 4*nWin/3 || std(x) == 0, continue; end
    x = x - mean(x);

    % (1) spectral concentration inside the band
    [P, F] = pwelch(x, hann(nWin), round(nWin/2), [], fs);
    P = P(:);  F = F(:);
    inB = F >= band(1) & F <= band(2);
    if ~any(inB) || sum(P(inB)) <= 0, continue; end
    Pm = P;  Pm(~inB) = -inf;
    [~, ip] = max(Pm);   f0(n) = F(ip);
    nearPk  = F >= f0(n)*0.8 & F <= f0(n)*1.2;
    conc(n) = sum(P(nearPk)) / sum(P(inB));

    % (2) autocorrelation side lobe of the band-passed trace (drift removed)
    xb = get_bandpass(x, fs, band);
    xb = xb(isfinite(xb));
    if numel(xb) < 3*maxLg, acPk(n) = 0; else
        ac = xcorr(xb - mean(xb), maxLg, 'coeff');
        ac = ac(maxLg+1:end);
        iz = find(ac <= 0, 1);
        if isempty(iz), iz = 2; end
        acPk(n) = max(ac(iz:end));
    end
    score(n) = max(conc(n), 0) * max(acPk(n), 0);
end
detail = struct('conc', conc, 'acPk', acPk);
end
