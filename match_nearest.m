function [idx_B, dist] = match_nearest(A, B)
% match_nearest  For each element of A, find the closest value in B.
%
% Usage:
%   [idx_B, dist] = match_nearest(A, B)
%
% Inputs:
%   A  - 1D vector of query values   (1 x Na  or  Na x 1)
%   B  - 1D vector of reference values (1 x Nb  or  Nb x 1)
%
% Outputs:
%   idx_B - Na x 1  index into B of the closest match for each element of A
%   dist  - Na x 1  signed distance  (A(i) - B(idx_B(i)))
%
% On an exact tie (A equidistant from two reference values) the larger
% reference value is chosen. NaN entries of A return NaN.
%
% Example:
%   A = [1.2, 3.7, 5.0];
%   B = [1.0, 2.5, 4.0, 6.0];
%   [idx, d] = match_nearest(A, B)
%   % idx = [1; 3; 4]   (B(1)=1.0, B(3)=4.0, B(4)=6.0 by tie rule)
%   % d   = [0.2; -0.3; -1.0]

A = A(:);   % ensure column
B = B(:);

% Nearest-neighbour via sorted B + midpoint bin edges (O(N log N) time, O(N) memory).
% NB: avoids the Na x Nb distance matrix, which is huge for long time axes.
if numel(B) == 1
    idx_B = ones(numel(A),1);
else
    [Bs, ord] = sort(B);
    edges = [-inf; (Bs(1:end-1)+Bs(2:end))/2; inf];   % A in [edge_k,edge_{k+1}) -> nearest Bs(k)
    k = discretize(A, edges);                          % NaN for NaN entries of A
    idx_B = nan(numel(A),1);
    good = ~isnan(k);
    idx_B(good) = ord(k(good));                        % map back to original (unsorted) B index
end

% Signed distance: positive means A is above the matched B value
dist = A - B(idx_B);

end
