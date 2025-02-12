function [avg_matrix centerEdge]= Bin_Vector(A, T, edges)
% SUM_MATRIX_BY_EDGES Sums matrix columns based on edge indices on the time axis.
%   summed_matrix = sum_matrix_by_edges(A, T, edges)
% 
% Parameters:
%   A      - Input matrix of size N x T.
%   T      - Time axis corresponding to the columns of A (1 x T).
%   edges  - Vector defining the time edges for summing (e.g., [0, 100, 200, ...]).
% 
% Returns:
%   summed_matrix - Matrix of size N x (length(edges)-1), where each column
%                   contains the summed values from the corresponding edge range.

% Validate inputs
if size(T, 1) > 1
    T = T'; % Ensure T is a row vector
end

if size(edges, 1) > 1
    edges = edges'; % Ensure edges is a row vector
end

if length(edges) < 2
    error('Edge vector must have at least two elements.');
end

% Ensure edges lie within the range of T
edges = max(edges, min(T));
edges = min(edges, max(T));

centerEdge=mean([edges(1:end-1); edges(2:end)],1);

% Find indices corresponding to edges on the T axis
[~, edge_indices] = histc(edges, [-Inf, T, Inf]);
edge_indices = max(edge_indices - 1, 1); % Ensure indices are valid

% Initialize output matrix
N = size(A, 1);
num_bins = length(edges) - 1;
avg_matrix = zeros(N, num_bins);

% Sum the columns for each edge range
for i = 1:num_bins
    col_start = edge_indices(i);
    col_end = edge_indices(i+1) - 1;
    col_end = min(col_end, size(A, 2)); % Ensure indices do not exceed bounds
    avg_matrix(:, i) = mean(A(:, col_start:col_end), 2,'omitnan');
end

end