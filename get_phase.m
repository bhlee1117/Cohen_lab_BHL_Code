function [cleanPhase, filteredTrace, Magnitude] = get_phase(trace, fs, Band, excludeMask)
% GET_PHASE  Instantaneous phase, band-limited trace, and envelope of a signal.
%
%   [phase, filtered, mag] = get_phase(trace, fs, Band)
%   Band-passes TRACE in the frequency range Band = [lo hi] (Hz), then uses the
%   Hilbert transform to return:
%       cleanPhase     instantaneous phase (rad, -pi..pi)      [column vector]
%       filteredTrace  band-passed trace                       [row vector]
%       Magnitude      analytic-signal amplitude (envelope)    [row vector]
%
%   [...] = get_phase(trace, fs, Band, excludeMask)
%   Treats the samples where EXCLUDEMASK is true (and any NaN samples) as
%   missing: they are interpolated away BEFORE filtering so that spike /
%   stimulus transients do not leak into the pass band, and the returned phase
%   and magnitude are set to NaN at those samples. Use this to keep post-spike
%   depolarizations from contaminating the phase estimate.
%
%   NOTE: filtfilt + hilbert are ACAUSAL -- the phase at any time point depends
%   on the whole TRACE, including samples AFTER that point. To estimate the
%   phase at a spike onset from PRE-spike dynamics only, use GET_PRESPIKE_PHASE,
%   which applies this function inside a causal pre-onset window.

    % Ensure trace is a column vector
    trace = trace(:);
    nT = numel(trace);

    if nargin < 3 || isempty(Band),        Band = [5 12]; end          % default theta band
    if nargin < 4 || isempty(excludeMask), excludeMask = false(nT,1); end
    excludeMask = excludeMask(:) | ~isfinite(trace);

    % Interpolate masked / NaN samples so that filtfilt (which cannot handle
    % NaN) and the Hilbert transform see a continuous, transient-free trace.
    x = trace;
    if any(excludeMask)
        good = ~excludeMask;
        if nnz(good) < 3
            cleanPhase = nan(nT,1); filteredTrace = nan(1,nT); Magnitude = nan(1,nT);
            return;
        end
        t = (1:nT)';
        x = interp1(t(good), trace(good), t, 'linear');
        x = fillmissing(x, 'nearest');   % fill any leading/trailing gaps
    end

    % Design a 2nd-order Butterworth band-pass and apply it zero-phase.
    [b, a] = butter(2, Band / (fs / 2), 'bandpass');
    filteredTrace = filtfilt(b, a, x);

    % Analytic signal via the Hilbert transform.
    analyticSignal = hilbert(filteredTrace);
    cleanPhase = angle(analyticSignal);
    Magnitude  = abs(analyticSignal);

    % Do not report phase / magnitude where the signal was actually missing.
    cleanPhase(excludeMask) = NaN;
    Magnitude(excludeMask)  = NaN;

    filteredTrace = filteredTrace(:)';
    Magnitude     = Magnitude(:)';
end
