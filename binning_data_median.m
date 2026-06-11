function [bin_medians, bin_stds, bin_centers, bin_membership, bin_values] = binning_data_median(data_cells, bin_edges)
% binning_data_median  Bin multi-recording data along an x-axis and compute median/std per bin.
%
% Identical to binning_data, but uses median instead of mean as the central
% tendency measure. Useful when per-bin distributions are skewed or contain
% outliers.
%
% Usage:
%   [bin_medians, bin_stds, bin_centers, bin_membership, bin_values] = ...
%       binning_data_median(data_cells, bin_edges)
%
% Inputs:
%   data_cells  - {nRecordings x 1} cell array. Each cell contains an
%                 [nSamples x 2] matrix where:
%                   column 1 = x-axis values (variable to bin by)
%                   column 2 = y-axis values (measure to summarize)
%                 Empty cells are skipped.
%   bin_edges   - [1 x (nBins+1)] vector of bin edges (left-inclusive,
%                 right-exclusive: bin j spans [edge(j), edge(j+1)) )
%
% Outputs:
%   bin_medians    - [1 x nBins] median of y-values in each bin (NaN if empty)
%   bin_stds       - [1 x nBins] std    of y-values in each bin (NaN if empty)
%   bin_centers    - [1 x nBins] center of each bin = mean(left, right edge)
%   bin_membership - {nBins x nRecordings} logical index vectors indicating
%                    which samples from each recording fall in each bin
%   bin_values     - {1 x nBins} cell array of all pooled y-values per bin
%
% Example:
%   bin_edges = 0:10:100;
%   [med, sd, xc] = binning_data_median(data_cells, bin_edges);
%   errorbar(xc, med, sd);
%
% See also: binning_data, binning_data2D
%
% 2024.10.27  Byung Hun Lee

nBins       = length(bin_edges) - 1;
nRecordings = length(data_cells);

bin_centers    = mean([bin_edges(1:end-1); bin_edges(2:end)], 1);
bin_medians    = nan(1, nBins);
bin_stds       = nan(1, nBins);
bin_membership = cell(nBins, nRecordings);
bin_values     = cell(1, nBins);

for iBin = 1:nBins
    pooled_y = [];   % accumulate y-values from all recordings for this bin

    for iRec = 1:nRecordings
        if isempty(data_cells{iRec})
            bin_membership{iBin, iRec} = [];
            continue
        end

        x_vals = data_cells{iRec}(:, 1);   % x-axis values for this recording
        y_vals = data_cells{iRec}(:, 2);   % y-axis values for this recording

        in_bin = (x_vals >= bin_edges(iBin)) & (x_vals < bin_edges(iBin + 1));
        bin_membership{iBin, iRec} = in_bin;

        pooled_y = [pooled_y, y_vals(in_bin)'];
    end

    if ~isempty(pooled_y)
        bin_medians(iBin) = median(pooled_y, 'omitnan');
        bin_stds(iBin)    = std(pooled_y,    'omitnan');
        bin_values{iBin}  = pooled_y;
    end
end
end
