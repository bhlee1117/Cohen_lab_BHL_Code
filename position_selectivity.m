function [p_map, z_map] = position_selectivity(V, pos, lap, varargin)
%POSITION_SELECTIVITY_VOLTAGE Compute position selectivity in raw voltage
%
%   [p_map, z_map] = position_selectivity_voltage(V, pos, t, lap, ...)
%
%   Optional name-value pairs:
%       'nBins'      - Number of spatial bins (default = 150)
%       'nShuffles'  - Number of circular shuffles (default = 1000)
%       'Show'       - 'on' (default) or 'off' to toggle figures

% --- Input validation ---
if nargin < 4
    error('Function requires at least four inputs: V, pos, t, and lap.');
end

if ~isvector(V) || ~isvector(pos) || ~isvector(lap)
    error('Inputs V, pos, t, and lap must be vectors.');
end

if length(V) ~= length(pos) || length(V) ~= length(lap)
    error('Inputs V, pos, t, and lap must have the same length.');
end

% --- Parse optional arguments ---
p = inputParser;
addParameter(p, 'nBins', 150, @(x) isnumeric(x) && isscalar(x) && x > 1);
addParameter(p, 'nShuffles', 1000, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, 'Show', 'on', @(x) ischar(x) || isstring(x));
parse(p, varargin{:});

nBins = p.Results.nBins;
nShuffles = p.Results.nShuffles;
showFigures = strcmpi(p.Results.Show, 'on');

% --- Interpolation onto uniform timebase ---
% Fs = 1000; % Hz
% t_sec = t / Fs;
% t_cont = (min(t):max(t)) / Fs;

validFrm = sum(isnan(V),1)==0;
V_int    = V(validFrm); %interp1(t_sec(validFrm), V(validFrm), t_cont, 'linear','extrap');
pos_int  = pos(validFrm); %interp1(t_sec, pos, t_cont, 'nearest', 'extrap');
lap_int  = lap(validFrm); %interp1(t_sec, lap, t_cont, 'nearest', 'extrap');

% --- Spatial binning ---
posEdges = linspace(min(pos_int), max(pos_int), nBins + 1);
posBins = discretize(pos_int, posEdges);

% --- Real map ---
voltageMap_real = nan(1, nBins);
for i = 1:nBins
    bin_idx = (posBins == i);
    if any(bin_idx)
        voltageMap_real(i) = mean(V_int(bin_idx), 'omitnan');  % Change to var, max, etc. if needed
    end
end

% --- Shuffling ---
voltageMap_shuffled = nan(nBins, nShuffles);
wb = waitbar(0, 'Shuffling voltage trace...');

for s = 1:nShuffles
    V_shuffled = nan(size(V_int));
    laps_unique = unique(lap_int);

    for l = 1:length(laps_unique)
        lap_id = laps_unique(l);
        lap_idx = find(lap_int == lap_id);
        if length(lap_idx) > 1
            shift_amt = randi(length(lap_idx));
            V_shuffled(lap_idx) = circshift(V_int(lap_idx), shift_amt);
        end
    end

    for i = 1:nBins
        bin_idx = (posBins == i);
        if any(bin_idx)
            voltageMap_shuffled(i, s) = mean(V_shuffled(bin_idx), 'omitnan');
        end
    end

    if mod(s, 10) == 0 || s == nShuffles
        waitbar(s / nShuffles, wb);
    end
end
close(wb);

% --- Z-score ---
shuff_mean = mean(voltageMap_shuffled, 2, 'omitnan');
shuff_std  = std(voltageMap_shuffled, 0, 2, 'omitnan');
z_map = (voltageMap_real' - shuff_mean) ./ shuff_std;

% --- P-value ---
p_map = mean(bsxfun(@ge, voltageMap_shuffled, voltageMap_real'), 2);
p_map(isnan(voltageMap_real')) = NaN;

% --- Optional plot ---
if showFigures
    figure;
    subplot(2,1,1);
    plot(1:nBins, z_map, 'k'); xlabel('Position Bin'); ylabel('Z-score');
    title('Z-scored Voltage Across Position');
    xlim([1 nBins]);

    subplot(2,1,2);
    plot(1:nBins, p_map, 'r'); xlabel('Position Bin'); ylabel('p-value');
    title('P-value from Shuffled Null Distribution');
    xlim([1 nBins]);
end
end
