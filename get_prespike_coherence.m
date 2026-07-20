function [pmi, dphi, ok] = get_prespike_coherence(traceA, traceB, fs, Band, onsets, excludeMask, W_ms, guard_ms)
% GET_PRESPIKE_COHERENCE  Causal pre-onset theta PHASE-MATCH index between two
% signals, e.g. distal vs basal subthreshold, using ONLY pre-onset samples so
% post-spike dynamics cannot bias the estimate.
%
%   [pmi, dphi, ok] = get_prespike_coherence(traceA, traceB, fs, Band, onsets, ...
%                                            excludeMask, W_ms, guard_ms)
%
% For each onset a pre-window (W_ms ending strictly before the onset) is taken
% from both traces, masked / NaN samples are interpolated away, both are
% band-passed (Band, Hz) and Hilbert-transformed, and the SIGNED phase-match
% index
%     pmi  = mean_t cos( phiA(t) - phiB(t) )
% is computed over the reliable interior of the window (both Hilbert edges
% dropped by guard_ms). pmi is +1 when A and B are IN PHASE (matched), 0 when
% orthogonal / incoherent, and -1 when ANTI-phase. dphi is the phaseA-phaseB
% difference read at (onset - guard).
%
% Defaults: W_ms = 600, guard_ms = 40. See also GET_PRESPIKE_PHASE, GET_PHASE.

    traceA = traceA(:)';  traceB = traceB(:)';
    nT = numel(traceA);
    if nargin < 4 || isempty(Band),        Band = [5 12]; end
    if nargin < 6 || isempty(excludeMask), excludeMask = false(1,nT); end
    if nargin < 7 || isempty(W_ms),        W_ms = 600; end
    if nargin < 8 || isempty(guard_ms),    guard_ms = 40; end

    excludeMask = excludeMask(:)' | ~isfinite(traceA) | ~isfinite(traceB);
    onsets = round(onsets(:)');
    K = numel(onsets);

    pmi  = nan(1,K);
    dphi = nan(1,K);
    ok   = false(1,K);

    W = round(W_ms   * fs/1000);
    g = round(guard_ms * fs/1000);
    minReliable = round(fs / min(Band));   % need >= ~1 cycle inside the window

    for k = 1:K
        t0 = onsets(k);
        s0 = t0 - W;
        if s0 < 1 || t0 > nT, continue; end

        idx = s0:(t0-1);
        m   = excludeMask(idx);
        if mean(m) > 0.5, continue; end
        good = ~m;
        if nnz(good) < 3, continue; end

        tt = 1:numel(idx);
        a = interp1(tt(good), traceA(idx(good)), tt, 'linear');  a = fillmissing(a,'nearest');
        b = interp1(tt(good), traceB(idx(good)), tt, 'linear');  b = fillmissing(b,'nearest');

        pa = get_phase(a - mean(a,'omitnan'), fs, Band);   % column phase
        pb = get_phase(b - mean(b,'omitnan'), fs, Band);

        rr = (g+1):(numel(idx)-g);      % reliable interior (drop both edges)
        if numel(rr) < minReliable, continue; end

        d = pa(rr) - pb(rr);
        pmi(k) = mean(cos(d), 'omitnan');   % signed: +1 in-phase, 0 orthogonal, -1 anti-phase

        r = numel(idx) - g;             % read point for the instantaneous offset
        dd = pa(r) - pb(r);
        dphi(k) = mod(dd + pi, 2*pi) - pi;
        ok(k) = true;
    end
end
