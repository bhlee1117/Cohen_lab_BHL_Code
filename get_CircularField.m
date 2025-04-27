function [fieldSizes, fieldEdges, fieldMat] = get_CircularField(sigMat, interval)
% get_CircularField: detects circular fields with merging and reconstructs index map
%
% Inputs:
%   sigMat   - N x T logical matrix (true = significance)
%   interval - max gap (bins) to merge fields (default = 0)
%
% Outputs:
%   fieldSizes  - N x 1 cell array: vector of merged field lengths
%   fieldEdges  - N x 1 cell array: [start, end] index pairs (1-based)
%   fieldMat    - N x T matrix: integer map labeling merged fields (0 = no field)

if nargin < 2
    interval = 0;
end

[N, T] = size(sigMat);
fieldSizes = cell(N, 1);
fieldEdges = cell(N, 1);
fieldMat = zeros(N, T);

for i = 1:N
    row = sigMat(i,:);
    labeled = bwlabel(row);
    stats = regionprops(labeled, 'PixelIdxList');
    bounds = cell2mat(cellfun(@(x) [x(1), x(end)], {stats.PixelIdxList}', 'UniformOutput', false));

    if isempty(bounds)
        fieldSizes{i} = [];
        fieldEdges{i} = [];
        continue;
    end

    % Step 1: sort and prepare bounds
    bounds = sortrows(bounds, 1);
    merged = bounds(1, :);

    % Step 2: merge based on interval
    for j = 2:size(bounds, 1)
        prevEnd = merged(end, 2);
        currStart = bounds(j, 1);
        gap = mod(currStart - prevEnd - 1, T);

        if gap <= interval
            merged(end, 2) = bounds(j, 2);
        else
            merged = [merged; bounds(j,:)]; %#ok<AGROW>
        end
    end

    % Step 3: wrap-around merge (last and first)
    if size(merged, 1) > 1
        gap = mod(merged(1,1) - merged(end,2) - 1, T);
        if gap <= interval
            merged(1,1) = merged(end,1);
            merged(end,:) = [];
        end
    end

    % Step 4: build outputs
    startIdx = mod(merged(:,1) - 1, T) + 1;
    endIdx   = mod(merged(:,2) - 1, T) + 1;
    lengths  = mod(endIdx - startIdx + T, T) + 1;

    fieldSizes{i} = lengths(:);
    fieldEdges{i} = [startIdx(:), endIdx(:)];

    % Step 5: reconstruct fieldMat with merged field indices
    for f = 1:size(merged,1)
        s = startIdx(f);
        e = endIdx(f);
        if s <= e
            fieldMat(i, s:e) = f;
        else
            fieldMat(i, [s:T, 1:e]) = f;
        end
    end
end
end
