function [SI_obs, SI_null] = SpatialInformation_shuffling(V_int, pos_int, lap_int, n_shuffles, nBins)
% Computes Skaggs' spatial information with lap-wise shuffling for null distribution.
%
% Inputs:
%   V_int      - Binary spike vector (1 x T)
%   pos_int    - Position vector (1 x T)
%   lap_int    - Lap ID vector (1 x T)
%   n_shuffles - Number of shuffles
%   nBins      - Number of spatial bins
%
% Outputs:
%   SI_obs     - Observed spatial information (bits/spike)
%   SI_null    - Null distribution

% --- Preprocess ---
valid = ~isnan(V_int) & ~isnan(pos_int) & ~isnan(lap_int);
V_int = V_int(valid);
pos_int = pos_int(valid);
lap_int = lap_int(valid);

% --- Compute observed SI ---
SI_obs = compute_SkaggsInfo(V_int, pos_int, nBins);

% --- Generate null distribution ---
SI_null = zeros(n_shuffles, 1);
laps = unique(lap_int);

h = waitbar(0, 'Lap-wise shuffling for Skaggs Info...');
for s = 1:n_shuffles
    V_shuffle = nan(size(V_int));
    for l = 1:length(laps)
        idx = find(lap_int == laps(l));
        if length(idx) > 1
            shift_amt = randi(length(idx));
            V_shuffle(idx) = circshift(V_int(idx), shift_amt);
        else
            V_shuffle(idx) = V_int(idx);
        end
    end

    SI_null(s) = compute_SkaggsInfo(V_shuffle, pos_int, nBins);
    SI_null(s) = max(SI_null(s), 0);  % ensure non-negative

    if mod(s, round(n_shuffles / 100)) == 0
        waitbar(s / n_shuffles, h);
    end
end
close(h);

fprintf('Observed Skaggs Info: %.4f bits/spike\n', SI_obs);
fprintf('p-value: %.4f\n', mean(SI_null >= SI_obs));


% ---------- Subfunction ----------
function SI = compute_SkaggsInfo(spikeVec, posVec, nBins)
    edges = linspace(min(posVec), max(posVec), nBins+1);
    bin_idx = discretize(posVec, edges);
    
    % Filter out invalid bins (e.g., outside edges)
    valid = ~isnan(bin_idx);
    bin_idx = bin_idx(valid);
    spikeVec = spikeVec(valid);
    
    occupancy = accumarray(bin_idx(:), 1, [nBins, 1]);
    spike_counts = accumarray(bin_idx(:), spikeVec(:), [nBins, 1]);

    p_i = occupancy / sum(occupancy);
    lambda_i = spike_counts ./ occupancy;
    lambda_i(isnan(lambda_i)) = 0;

    lambda = sum(spike_counts) / sum(occupancy);

    info_per_bin = p_i .* (lambda_i / lambda) .* log2(lambda_i / lambda);
    info_per_bin(~isfinite(info_per_bin)) = 0;

    SI = sum(info_per_bin);
end


end
