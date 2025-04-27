function [p_map, z_map, freqEdges] = position_frequency_selectivity(V, pos, t, lap, varargin)
%POSITION_SELECTIVITY Compute and assess position selectivity in voltage power
%
%   [p_map, z_map, freqEdges] = position_frequency_selectivity(V, pos, t, lap, ...)
%
%   Required Inputs:
%       V        - 1xN voltage trace
%       pos      - 1xN position vector
%       t        - 1xN time vector (in ms)
%       lap      - 1xN lap index vector
%
%   Optional name-value pairs:
%       'nBins'      - Number of spatial bins (default = 150)
%       'nShuffles'  - Number of circular shuffles (default = 1000)
%       'freqEdges'  - Frequency bin edges (default = 0:4:100)
%       'Show'       - 'on' (default) or 'off' to toggle figures
%
%   Example usage:
%       position_selectivity(V, pos, t, lap, 'nBins', 100, 'Show', 'off');

if nargin < 4
    error('Function requires at least four inputs: V, pos, t, and lap.');
end

% Check vector and size consistency
if ~isvector(V) || ~isvector(pos) || ~isvector(t) || ~isvector(lap)
    error('Inputs V, pos, t, and lap must be vectors.');
end

if length(V) ~= length(pos) || length(V) ~= length(t) || length(V) ~= length(lap)
    error('Inputs V, pos, t, and lap must have the same length.');
end

% --- Parse optional arguments ---
p = inputParser;
addParameter(p, 'nBins', 150, @(x) isnumeric(x) && isscalar(x) && x > 1);
addParameter(p, 'nShuffles', 1000, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, 'freqEdges', 0:4:100, @(x) isnumeric(x) && isvector(x) && all(diff(x) > 0));
addParameter(p, 'Show', 'on', @(x) ischar(x) || isstring(x));
parse(p, varargin{:});

% Assign to variables
nBins      = p.Results.nBins;
nShuffles  = p.Results.nShuffles;
freqEdges  = p.Results.freqEdges;
showFigures = strcmpi(p.Results.Show, 'on');

% --- Parameters ---
Fs = 1000;  % Hz
window_length = 512;
noverlap = 500;
nfft = 256;

% --- Interpolation onto uniform time ---
t_sec = t / Fs; %in original axis
t_cont = (min(t):max(t)) / Fs;
validFrm = ~isnan(V);
validT = t(validFrm);

V_int    = interp1(t_sec(validFrm), V(validFrm), t_cont, 'linear','extrap');
pos_int  = interp1(t_sec, pos, t_cont, 'nearest', 'extrap');
lap_int  = interp1(t_sec, lap, t_cont, 'nearest', 'extrap');

% --- Spectrogram on original data ---
[S, F, T_spec] = spectrogram(V_int, window_length, noverlap, nfft, Fs);
T_spec=T_spec+t_cont(1);
P = abs(S).^2;
nFreqBins = length(freqEdges) - 1;
freqBinIdx = discretize(F, freqEdges);  % assign each frequency to a bin

% --- Map time points to position bins ---
pos_T = interp1(t_cont, pos_int, T_spec, 'nearest', 'extrap');
lap_T = interp1(t_cont, lap_int, T_spec, 'nearest', 'extrap');
posEdges = linspace(min(pos_T), max(pos_T), nBins + 1);
posBin_T = discretize(pos_T, posEdges);
posBin_T = posBin_T(ismember(T_spec*Fs,validT));

lap_T=lap_T(:,ismember(T_spec*Fs,validT));
P=P(:,ismember(T_spec*Fs,validT));

nFreqs = length(F);
powerMap_real = nan(nFreqBins, nBins);

for i = 1:nBins
    bin_idx = (posBin_T == i);
    if any(bin_idx)
        for f = 1:nFreqBins
            freq_idx = (freqBinIdx == f);
            if any(freq_idx)
                powerMap_real(f, i) = mean(P(freq_idx, bin_idx), 'all', 'omitnan');
            end
        end
    end
end

% --- Shuffling ---
powerMap_shuffled = nan(nFreqBins, nBins, nShuffles);
wb = waitbar(0, 'Shuffling position labels...');

for s = 1:nShuffles
    shift_amt = randi(length(posBin_T));
    posBin_shuffled = nan(size(posBin_T));

    laps_unique = unique(lap_T);
    for l = 1:length(laps_unique)
        lap_id = laps_unique(l);
        lap_idx = find(lap_T == lap_id);

        if length(lap_idx) > 1
            shift_amt = randi(length(lap_idx));
            posBin_shuffled(lap_idx) = circshift(posBin_T(lap_idx), shift_amt);
        end
    end


    for i = 1:nBins
        bin_idx = (posBin_shuffled == i);
        if any(bin_idx)
            for f = 1:nFreqBins
                freq_idx = (freqBinIdx == f);
                if any(freq_idx)
                    powerMap_shuffled(f, i, s) = mean(P(freq_idx, bin_idx), 'all', 'omitnan');
                end
            end
        end
    end
    if mod(s,10)==0 || s == nShuffles
        waitbar(s / nShuffles, wb);
    end
end
close(wb);


% --- Z-score ---
shuff_mean = mean(powerMap_shuffled, 3, 'omitnan');
shuff_std  = std(powerMap_shuffled, 0, 3, 'omitnan');
z_map = (powerMap_real - shuff_mean) ./ shuff_std;

% --- P-value ---
p_map = mean(bsxfun(@ge, powerMap_shuffled, powerMap_real), 3);  % observed greater than shuffle
p_map(isnan(powerMap_real)) = NaN;  % mask empty bins

% --- Plot Z-score Map ---
if showFigures
figure; tiledlayout(2,1);
ax2=nexttile([1 1]);
imagesc(1:nBins, freqEdges(1:end-1) + diff(freqEdges)/2, z_map);
axis xy;
xlabel('Position Bin');
ylabel('Frequency (Hz)');
title('Z-scored Power Across Position');
colormap turbo; colorbar;

% --- Optional: Thresholded map ---
ax1=nexttile([1 1]);
imagesc(1:nBins, freqEdges(1:end-1) + diff(freqEdges)/2, p_map);
axis xy;
xlabel('Position Bin');
ylabel('Frequency (Hz)');
title('Significant Power (p < 0.01)');
colormap(ax1,flipud(turbo(256)))
end
end