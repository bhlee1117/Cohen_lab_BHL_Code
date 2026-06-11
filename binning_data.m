function result = binning_data(data_cells, bin_edges)
% binning_data  Bin multi-recording data along an x-axis and compute statistics per bin.
%
% For each bin defined by bin_edges, collects all y-values across every
% recording whose x-value falls in that bin, then returns a struct of
% summary statistics.
%
% Usage:
%   result = binning_data(data_cells, bin_edges)
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
% Output:
%   result  - scalar struct with fields:
%     .centers     [1 x nBins]          bin center positions
%     .mean        [1 x nBins]          mean   of pooled y-values (NaN if empty)
%     .std         [1 x nBins]          std    of pooled y-values (NaN if empty)
%     .sem         [1 x nBins]          SEM    of pooled y-values (NaN if empty)
%                                         = std / sqrt(nRecordings with data in bin)
%     .values      {1 x nBins}          pooled y-values per bin
%     .membership  {nBins x nRecordings} logical vectors: which samples per
%                                         recording fall in each bin
%
% Example:
%   bin_edges = 0:10:100;
%   result = binning_data(data_cells, bin_edges);
%   errorbar(result.centers, result.mean, result.sem);
%
% See also: binning_data_median, binning_data2D
%
% 2024.10.27  Byung Hun Lee

nBins       = length(bin_edges) - 1;
nRecordings = length(data_cells);

% Preallocate output struct
result.centers    = mean([bin_edges(1:end-1); bin_edges(2:end)], 1);
result.mean       = nan(1, nBins);
result.std        = nan(1, nBins);
result.sem        = nan(1, nBins);
result.values     = cell(1, nBins);
result.membership = cell(nBins, nRecordings);

for iBin = 1:nBins
    pooled_y = [];   % accumulate y-values from all recordings for this bin

    for iRec = 1:nRecordings
        if isempty(data_cells{iRec})
            result.membership{iBin, iRec} = [];
            continue
        end

        x_vals = data_cells{iRec}(:, 1);   % x-axis values for this recording
        y_vals = data_cells{iRec}(:, 2);   % y-axis values for this recording

        in_bin = (x_vals >= bin_edges(iBin)) & (x_vals < bin_edges(iBin + 1));
        result.membership{iBin, iRec} = in_bin;

        pooled_y = [pooled_y, y_vals(in_bin)'];
    end

    if ~isempty(pooled_y)
        % nRec_in_bin: number of recordings that contributed data to this bin
        nRec_in_bin = sum(cellfun(@any, result.membership(iBin, :)));

        result.mean(iBin)   = mean(pooled_y, 'omitnan');
        result.std(iBin)    = std(pooled_y,  'omitnan');
        result.sem(iBin)    = result.std(iBin) / sqrt(nRec_in_bin);
        result.values{iBin} = pooled_y;
    end
end
end
