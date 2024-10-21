function spatial_info = SpatialInfo(firing_rate_maps, occupancy)

   % Sum the firing rate maps and occupancy times across all laps
    total_spike_counts = nansum(firing_rate_maps .* occupancy, 1);  % Total spike counts per bin
    total_occupancy = nansum(occupancy, 1);  % Total occupancy per bin

    % Replace zero occupancy with NaN to avoid division by zero
    total_occupancy(total_occupancy == 0) = NaN;

    % Calculate the probability of visiting each bin (normalized occupancy)
    p_i = total_occupancy / nansum(total_occupancy);

    % Calculate the mean firing rate per bin
    r_i = total_spike_counts ./ total_occupancy;

    % Calculate the overall mean firing rate across all bins
    r = nansum(r_i .* p_i);  % Weighted average firing rate

    % Compute the spatial information for the pooled data
    spatial_info = sum(p_i .* (r_i / r) .* log2(r_i / r),2, 'omitnan');
end