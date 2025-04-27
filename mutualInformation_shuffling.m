function [MI_obs, MI_null] = mutualInformation_shuffling(V_int, pos_int, lap_int, n_shuffles)
% Computes mutual information between interpolated voltage and position
% using lap-wise circular shuffling.
% Inputs:
%   V_int      - Interpolated voltage (1 x T)
%   pos_int    - Interpolated position (1 x T)
%   lap_int    - Interpolated lap IDs (1 x T)
%   n_shuffles - Number of permutations
%
% Output:
%   MI_obs     - Observed mutual information (in bits)
%   MI_null    - Null distribution of MI under shuffled condition

L = 256;  % Number of histogram bins

% Remove NaNs
valid_idx = ~isnan(V_int) & ~isnan(pos_int) & ~isnan(lap_int);
V_int = V_int(valid_idx);
pos_int = pos_int(valid_idx);
lap_int = lap_int(valid_idx);

% Observed MI
MI_obs = mi(V_int(:), pos_int(:), L);
%MI_obs = mi_cont_cont(V_int(:), pos_int(:));

% Shuffle
MI_null = zeros(n_shuffles, 1);
laps_unique = unique(lap_int);

h = waitbar(0, 'Lap-wise shuffling for MI...');

for s = 1:n_shuffles
    V_shuffled = nan(size(V_int));

    for l = 1:length(laps_unique)
        idx = find(lap_int == laps_unique(l));
        if length(idx) > 1
            shift_amt = randi(length(idx));
            V_shuffled(idx) = circshift(V_int(idx), shift_amt);
        else
            V_shuffled(idx) = V_int(idx);
        end
    end

    MI_null(s) = mi(V_shuffled(:), pos_int(:), L);

    if mod(s, round(n_shuffles / 100)) == 0
        waitbar(s / n_shuffles, h);
    end
end
close(h);

fprintf('Observed MI: %.4f bits\n', MI_obs);
fprintf('p-value: %.4f\n', mean(MI_null >= MI_obs));
end
