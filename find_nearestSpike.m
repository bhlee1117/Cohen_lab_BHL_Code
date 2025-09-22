function [before_frame, after_frame, before_dist, after_dist] = find_nearestSpike(idx, spikevec)
% [before_frame, after_frame, before_dist, after_dist] = find_nearestSpike(idx, spikevec)
%   Finds the nearest 1 before and after each index in a binary vector.
%
% Inputs:
%   idx     - vector of indices
%   binvec  - binary vector (1D array)
%
% Outputs:
%   before_frame - nearest index <= idx(i) where binvec == 1
%   after_frame  - nearest index >= idx(i) where binvec == 1
%   before_dist  - idx(i) - before_frame
%   after_dist   - after_frame - idx(i)

% Preallocate
N = length(idx);
before_frame = NaN(size(idx));
after_frame  = NaN(size(idx));
before_dist  = Inf(size(idx));
after_dist   = Inf(size(idx));

% Find all 1s in binvec
one_locs = find(spikevec);

for i = 1:N
    current = idx(i);

    % Check if current index is already 1
    if spikevec(current) == 1
        before_frame(i) = current;
        after_frame(i)  = current;
        before_dist(i)  = 0;
        after_dist(i)   = 0;
        continue;
    end

    % Find nearest 1 before or at
    past_ones = one_locs(one_locs <= current);
    if ~isempty(past_ones)
        before_frame(i) = past_ones(end);
        before_dist(i) = current - before_frame(i);
    end

    % Find nearest 1 after or at
    future_ones = one_locs(one_locs >= current);
    if ~isempty(future_ones)
        after_frame(i) = future_ones(1);
        after_dist(i) = after_frame(i) - current;
    end
end
