function [binned_matrix, bin_centers] = bin_mat(dax, A, bin_edges)
% bin_mat Bins rows of matrix A based on dax values.
%
%   [binned_matrix, bin_centers] = bin_mat(dax, A, bin_edges)
%
%   Inputs:
%       dax        - 1 x N vector (binning axis)
%       A          - N x T matrix (data to be binned)
%       bin_edges  - 1 x R vector of bin edges
%
%   Outputs:
%       binned_matrix - (R-1) x T matrix, each row is the mean of rows in A within a bin
%       bin_centers   - 1 x (R-1) vector of bin centers

    % Validate input dimensions
    assert(isvector(dax) && isrow(dax), 'dax must be a 1 x N row vector');
    assert(size(A, 1) == length(dax), 'Number of rows in A must match length of dax');
    assert(isvector(bin_edges) && length(bin_edges) >= 2, 'bin_edges must have at least 2 elements');

    % Calculate bin centers
    bin_centers = mean([bin_edges(1:end-1); bin_edges(2:end)], 1);
    num_bins = length(bin_centers);
    T = size(A, 2);

    % Initialize output
    binned_matrix = NaN(num_bins, T);

    % Bin and average
    for i = 1:num_bins
        in_bin = dax >= bin_edges(i) & dax < bin_edges(i+1);
        if any(in_bin)
            binned_matrix(i, :) = mean(A(in_bin, :), 1, 'omitnan');
        end
    end
end
