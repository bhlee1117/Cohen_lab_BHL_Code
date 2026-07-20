function [phi, mag, finst, ok] = get_prespike_phase(trace, fs, Band, onsets, excludeMask, W_ms, guard_ms)
% GET_PRESPIKE_PHASE  Band-limited (theta) phase at event onsets, estimated
% CAUSALLY from PRE-onset data only, so that the post-spike depolarization can
% never bias the estimate. This answers "at which theta phase (measured from the
% pre-spike dynamics) was the spike triggered?".
%
%   [phi, mag, finst, ok] = get_prespike_phase(trace, fs, Band, onsets, ...
%                                               excludeMask, W_ms, guard_ms)
%
% INPUT
%   trace       1 x nTime signal. Use the spike-removed subthreshold trace.
%   fs          sampling rate (Hz)
%   Band        [lo hi] pass band (Hz), e.g. [5 12] for theta
%   onsets      1 x K frame indices (spike / event onsets) to report phase at
%   excludeMask 1 x nTime logical, true = ignore this sample (earlier spikes,
%               blue stim, NaN). Masked samples are interpolated away inside
%               each pre-window so their transients do not leak into the band.
%               (Optional; default: only non-finite samples of TRACE.)
%   W_ms        pre-window length (ms), default 600 (several theta cycles)
%   guard_ms    gap (ms) just before the onset that is skipped to avoid the
%               Hilbert end-edge artefact, default 40. The phase and the local
%               instantaneous frequency are read at (onset - guard) where the
%               estimate is reliable, then the phase is extrapolated forward
%               across the guard to the onset. Nothing after the onset is used.
%
% OUTPUT (all 1 x K, aligned to ONSETS)
%   phi     theta phase at the onset (rad, -pi..pi); NaN where not estimable
%   mag     theta envelope (amplitude) at onset-guard; NaN where not estimable
%   finst   local instantaneous theta frequency used for extrapolation (Hz)
%   ok      logical, true where a valid estimate was produced
%
% See also GET_PHASE.

    trace = trace(:)';
    nT = numel(trace);

    if nargin < 3 || isempty(Band),     Band = [5 12]; end
    if nargin < 5 || isempty(excludeMask), excludeMask = false(1,nT); end
    if nargin < 6 || isempty(W_ms),     W_ms = 600; end
    if nargin < 7 || isempty(guard_ms), guard_ms = 40; end

    excludeMask = excludeMask(:)' | ~isfinite(trace);
    onsets = round(onsets(:)');
    K = numel(onsets);

    phi   = nan(1,K);
    mag   = nan(1,K);
    finst = nan(1,K);
    ok    = false(1,K);

    W = round(W_ms   * fs/1000);   % pre-window length (samples)
    g = round(guard_ms * fs/1000); % edge guard (samples)
    Lf = round(fs / min(Band));    % ~one slow cycle, used for the freq estimate
    minReliable = round(fs / min(Band));   % need >= ~1 cycle before the read point

    for k = 1:K
        t0 = onsets(k);
        segStart = t0 - W;
        if segStart < 1 || t0 > nT, continue; end

        idx = segStart:(t0-1);         % strictly BEFORE the onset -> causal
        seg = trace(idx);
        m   = excludeMask(idx);
        if mean(m) > 0.5, continue; end % too much of the window is masked

        % interpolate masked / NaN samples inside the window
        good = ~m;
        if nnz(good) < 3, continue; end
        tt  = 1:numel(seg);
        seg = interp1(tt(good), seg(good), tt, 'linear');
        seg = fillmissing(seg, 'nearest');
        seg = seg - mean(seg, 'omitnan');

        % band-pass + analytic signal (filtfilt sees ONLY pre-onset samples)
        [ph, ~, mg] = get_phase(seg, fs, Band);   % ph column, mg row

        % reliable read point: guard samples before the window end
        r = numel(seg) - g;
        if r < minReliable, continue; end

        % local instantaneous frequency over ~one cycle ending at r
        rr = max(1, r-Lf):r;
        dphi = diff(unwrap(ph(rr)));
        fi = median(dphi, 'omitnan') * fs / (2*pi);
        if ~isfinite(fi) || fi <= 0, fi = mean(Band); end

        % extrapolate the phase from r forward to the onset (g+1 samples)
        nStep = g + 1;
        phRaw = ph(r) + 2*pi*fi/fs*nStep;
        phi(k)   = mod(phRaw + pi, 2*pi) - pi;    % wrap to (-pi, pi]
        mag(k)   = mg(r);
        finst(k) = fi;
        ok(k)    = true;
    end
end
