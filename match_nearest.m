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
% Example:
%   A = [1.2, 3.7, 5.0];
%   B = [1.0, 2.5, 4.0, 6.0];
%   [idx, d] = match_nearest(A, B)
%   % idx = [1; 3; 3]   (B(1)=1.0, B(3)=4.0, B(3)=4.0)
%   % d   = [0.2; -0.3; 1.0]

A = A(:);   % ensure column
B = B(:);

% Compute absolute distance matrix (Na x Nb) and find minimum along B axis
[~,        idx_B] = min(abs(A - B'), [], 2);   % Na x 1

% Signed distance: positive means A is above the matched B value
dist = A - B(idx_B);

end
