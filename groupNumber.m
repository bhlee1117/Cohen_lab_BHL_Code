function shrunkVec = groupNumber(vec, threshold)
% groupNumber Shrinks a vector by replacing clusters of values within a threshold with their mean.
%
%   shrunkVec = groupNumber(vec, threshold) sorts the input vector 'vec', groups consecutive
%   values closer than 'threshold', replaces each group with its mean, and returns the shrunk vector.
%
%   Input:
%       vec       - Input vector (row or column).
%       threshold - Minimum distance to keep values separate.
%
%   Output:
%       shrunkVec - Shrunk vector with means of close groups.
%
%   Notes: Handles empty vectors. Output matches input orientation.
if isempty(vec)
    shrunkVec = [];
    return;
end
sortedVec = sort(vec(:));  % Sort and ensure column vector
% Initialize groups
groups = {sortedVec(1)};
for i = 2:length(sortedVec)
    if abs(sortedVec(i) - groups{end}(end)) < threshold
        groups{end}(end+1) = sortedVec(i);  % Add to current group
    else
        groups{end+1} = sortedVec(i);  % Start new group
    end
end

% Compute means
shrunkVec = cellfun(@mean, groups);
if size(vec,1) == 1  % If input was row vector
shrunkVec = shrunkVec';
end
end