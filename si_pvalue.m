function [SI, SI_p, SI_null] = si_pvalue(spikeMat, binTr, tax, nShuffle, seed, label)
% SI_PVALUE  Spatial information + circular-shuffle significance.
%   [SI, SI_p, SI_null] = si_pvalue(spikeMat, binTr, tax, nShuffle, seed, label)
%   spikeMat : nROI x nT  spike train(s)
%   binTr    : 1 x nT  position-bin index per frame, NaN = excluded (non-running)
%   tax      : 1 x nT  time axis (s); used only to set a >= ~1 s minimum shift
%   nShuffle : # circular shuffles                 seed : rng seed (default 0)
%   label    : optional string; if given, prints shuffle progress + ETA
%
%   SI     : nROI x 1  observed spatial information (bits/spike)
%   SI_p   : nROI x 1  right-tailed, add-one-corrected p-value from the null
%   SI_null: nROI x nShuffle  null SI values
% Each shuffle rolls the spike train by a random time offset (binTr fixed), so
% spike count/autocorrelation are preserved while spike-position coupling is broken.
if nargin < 5 || isempty(seed), seed = 0; end
if nargin < 6, label = ''; end
verbose = ~isempty(label);
nT = size(spikeMat,2);
SI = spatialInfo_mat(spikeMat, binTr);

dt       = median(diff(tax(1:nT)));
minShift = max(1, round(1/dt));                       % >= ~1 s to break local autocorrelation
minShift = min(minShift, max(1,floor(nT/2)));
hiShift  = max(minShift+1, nT-minShift);

rng(seed);                                            % reproducible null
shifts   = randi([minShift, hiShift], 1, nShuffle);
SI_null  = zeros(size(spikeMat,1), nShuffle);
if verbose, t0 = tic; step = max(1, round(nShuffle/10)); end
for i = 1:nShuffle
    SI_null(:,i) = spatialInfo_mat(circshift(spikeMat, shifts(i), 2), binTr);
    if verbose && (mod(i,step)==0 || i==nShuffle)
        el = toc(t0);
        fprintf('    [%s] shuffle %4d/%d (%3.0f%%) | %.1fs elapsed, ~%.1fs left\n', ...
            label, i, nShuffle, 100*i/nShuffle, el, el*(nShuffle-i)/i);
    end
end
SI_p = (1 + sum(SI_null >= SI, 2)) ./ (nShuffle + 1);
SI_p(isnan(SI)) = NaN;                                % undefined SI (e.g. no spikes)
end
