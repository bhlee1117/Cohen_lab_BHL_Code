%function x_clean = remove_artifact_stft(x, fs, f_artifact, varargin)
function[x_clean, noiseFreqs, removed, info] = remove_artifact_stft(x, Fs, protectedBands,varargin)
% removePersistentPeriodicNoise
%
% Finds persistent periodic noise in a signal and removes it with notch filters,
% while avoiding user-defined protected frequency bands.
%
% Inputs
% ------
% x              : input signal, vector
% Fs             : sampling rate in Hz
% protectedBands : Nx2 matrix of protected frequency ranges in Hz
%                  Example: [0 2; 8 12] protects 0-2 Hz and 8-12 Hz
%
% Optional name-value inputs
% --------------------------
% 'FreqRange'        : frequency search range, default [0 Fs/2]
% 'WindowSec'        : spectrogram window size in seconds, default 2
% 'OverlapFrac'      : spectrogram overlap fraction, default 0.75
% 'PersistenceThr'   : fraction of time bins where noise must appear, default 0.6
% 'ProminenceFactor' : peak threshold strength, default 4
% 'Q'                : notch filter quality factor, default 20
% 'MinPeakDistanceHz': minimum distance between detected peaks, default 1
% 'PlotFlag'         : true/false, default true
%
% Outputs
% -------
% x_clean   : filtered signal
% noiseFreqs: detected noise frequencies in Hz
% removed   : x - x_clean
% info      : diagnostic structure

% -------------------------
% Parse inputs
% -------------------------
p = inputParser;
addParameter(p, 'FreqRange', [0 Fs/2]);
addParameter(p, 'WindowSec', 2);
addParameter(p, 'OverlapFrac', 0.75);
addParameter(p, 'PersistenceThr', 0.6);
addParameter(p, 'ProminenceFactor', 4);
addParameter(p, 'Q', 20);
addParameter(p, 'MinPeakDistanceHz', 1);
addParameter(p, 'PlotFlag', true);
parse(p, varargin{:});

freqRange        = p.Results.FreqRange;
windowSec        = p.Results.WindowSec;
overlapFrac      = p.Results.OverlapFrac;
persistenceThr   = p.Results.PersistenceThr;
prominenceFactor = p.Results.ProminenceFactor;
Q                = p.Results.Q;
minPeakDistHz    = p.Results.MinPeakDistanceHz;
plotFlag         = p.Results.PlotFlag;

x = x(:);
x = x - median(x, 'omitnan');

% Replace NaNs if present
if any(isnan(x))
    x = fillmissing(x, 'linear', 'EndValues', 'nearest');
end

% -------------------------
% Spectrogram
% -------------------------
winLength = round(windowSec * Fs);
winLength = max(winLength, 16);

noverlap = round(overlapFrac * winLength);
nfft = max(2^nextpow2(winLength), winLength);

[S, F, T] = spectrogram(x, hann(winLength), noverlap, nfft, Fs);

P = abs(S).^2;

% Keep only requested frequency range
freqMask = F >= freqRange(1) & F <= freqRange(2);
F2 = F(freqMask);
P2 = P(freqMask, :);

% -------------------------
% Exclude protected bands
% -------------------------
allowedMask = true(size(F2));

for k = 1:size(protectedBands, 1)
    allowedMask = allowedMask & ...
        ~(F2 >= protectedBands(k,1) & F2 <= protectedBands(k,2));
end

% Also avoid DC
allowedMask = allowedMask & F2 > 0;

% -------------------------
% Detect persistent narrowband peaks
% -------------------------
% Normalize each time bin to reduce non-stationary amplitude effects
Pnorm = P2 ./ median(P2, 1, 'omitnan');

% Robust frequency score over time
freqScore = median(Pnorm, 2, 'omitnan');

% Persistence: how often each frequency is unusually strong
localBaseline = median(Pnorm, 1, 'omitnan');
strongMask = Pnorm > prominenceFactor * localBaseline;
persistence = mean(strongMask, 2, 'omitnan');

% Only consider allowed frequencies
scoreForPeaks = freqScore;
scoreForPeaks(~allowedMask) = 0;
scoreForPeaks(persistence < persistenceThr) = 0;

% Peak threshold
noiseFloor = median(scoreForPeaks(scoreForPeaks > 0), 'omitnan');
if isempty(noiseFloor) || isnan(noiseFloor)
    noiseFreqs = [];
    x_clean = x;
    removed = zeros(size(x));
    info.F = F2;
    info.T = T;
    info.P = P2;
    info.freqScore = freqScore;
    info.persistence = persistence;
    return;
end

minPeakDistanceBins = max(1, round(minPeakDistHz / median(diff(F2))));

[~, locs] = findpeaks(scoreForPeaks, ...
    'MinPeakDistance', minPeakDistanceBins, ...
    'MinPeakProminence', noiseFloor);

noiseFreqs = F2(locs);

% Make sure detected frequencies are not inside protected bands
keep = true(size(noiseFreqs));
for k = 1:size(protectedBands, 1)
    keep = keep & ~(noiseFreqs >= protectedBands(k,1) & noiseFreqs <= protectedBands(k,2));
end
noiseFreqs = noiseFreqs(keep);

% -------------------------
% Apply notch filters
% -------------------------
x_clean = x;

for k = 1:numel(noiseFreqs)
    f0 = noiseFreqs(k);

    if f0 <= 0 || f0 >= Fs/2
        continue;
    end

    wo = f0 / (Fs/2);
    bw = wo / Q;

    [b, a] = iirnotch(wo, bw);

    x_clean = filtfilt(b, a, x_clean);
end

removed = x - x_clean;

% -------------------------
% Store diagnostics
% -------------------------
info.F = F2;
info.T = T;
info.P = P2;
info.freqScore = freqScore;
info.persistence = persistence;
info.allowedMask = allowedMask;
info.scoreForPeaks = scoreForPeaks;
info.noiseFreqs = noiseFreqs;
info.Q = Q;
info.protectedBands = protectedBands;

% -------------------------
% Plot diagnostics
% -------------------------
if plotFlag
    t = (0:numel(x)-1) / Fs;

    figure;

    subplot(4,1,1);
    plot(t, x);
    title('Raw signal');
    xlabel('Time (s)');
    ylabel('Amplitude');

    subplot(4,1,2);
    imagesc(T, F2, 10*log10(P2 + eps));
    axis xy;
    title('Spectrogram');
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    colorbar;

    subplot(4,1,3);
    plot(F2, freqScore);
    hold on;
    plot(F2, persistence * max(freqScore), '--');
    xline(noiseFreqs, 'r');
    title('Frequency score and detected persistent noise');
    xlabel('Frequency (Hz)');
    ylabel('Score');
    legend('Median normalized power', 'Persistence scaled', 'Detected noise');

    subplot(4,1,4);
    plot(t, x_clean);
    title('Filtered signal');
    xlabel('Time (s)');
    ylabel('Amplitude');
end

end
