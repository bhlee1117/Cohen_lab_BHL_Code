function [z_scored_data, mean_sampled_data, x_bin_centers, indicies] = zscore_binning(Data, x_bin_edges)
% [z_scored_data, mean_sampled_data, x_bin_centers, indicies] = zscore_binning(Data, x_bin_edges)
% - Data: cell array, each cell is N x 2 matrix [x, amplitude]
% - x_bin_edges: vector of bin edges
% - z_scored_data: cell array of [#bins x 2] matrices: [x_bin_center, zscore]

x_bin_centers = mean([x_bin_edges(1:end-1); x_bin_edges(2:end)],1);
num_bins = length(x_bin_centers);
num_samples = length(Data);

sample_means = NaN(num_samples, num_bins);  % Rows: samples, Cols: bins
indicies = cell(num_bins, num_samples);     % Store bin indices per sample and bin

% Step 1: compute sample-wise mean per bin using cellfun
for j = 1:num_bins
    bin_start = x_bin_edges(j);
    bin_end = x_bin_edges(j+1);

    % Apply binning logic with cellfun
    is_in_bin = cellfun(@(d) (d(:,1) >= bin_start) & (d(:,1) < bin_end), Data, 'UniformOutput', false);
    bin_vals = cellfun(@(d, idx) mean(d(idx,2), 'omitnan'), Data, is_in_bin, 'UniformOutput', true);

    sample_means(:, j) = bin_vals(:);
    indicies(j, :) = is_in_bin;
end

% Step 2: compute z-score across samples (per bin)
mu = mean(sample_means, 1, 'omitnan');
sigma = std(sample_means, 0, 1, 'omitnan');
z_scores = (sample_means - mu) ./ sigma;


% Step 3: package into cell array: each cell is [x_bin_centers, zscores]
z_scored_data = cell(num_samples, 1);
mean_sampled_data=cell(num_samples, 1);
for f = 1:num_samples
    z_scored_data{f} = [x_bin_centers(:), z_scores(f, :)'];
    mean_sampled_data{f} = [x_bin_centers(:), sample_means(f, :)'];
end

end
