function [roiOut, mapping] = get_combinedROI(roiIn, minSize)
% Combines (merges) small ROI masks in roiIn so that each output ROI has
% at least minSize pixels (whenever possible).
%
% Strategy:
%   1) Break each input slice into connected components (like your splitter).
%   2) Treat each component as an ROI candidate.
%   3) While there exists an ROI < minSize, merge it with its nearest neighbor
%      (by centroid distance) to form a larger ROI.
%
% Notes / assumptions:
%   - ROIs are assumed mostly non-overlapping. If overlaps occur, values are
%     combined using max() at overlapping pixels.
%   - If total pixels across ALL ROIs is < minSize, it is impossible to make
%     an ROI >= minSize; you will end up with a single ROI smaller than minSize.
%
% INPUTS:
%   roiIn   - X x Y x N matrix. Each slice is an ROI with non-zero values.
%   minSize - Minimum desired number of pixels for each output ROI.
%
% OUTPUTS:
%   roiOut  - X x Y x M matrix. Each slice is a merged ROI with preserved values.
%   mapping - M x 1 cell array. mapping{i} lists which original ROI slices (1..N)
%             contributed pixels to roiOut(:,:,i).

    if nargin < 2
        error('Usage: [roiOut, mapping] = get_combinedROI(roiIn, minSize)');
    end
    if minSize <= 0
        error('minSize must be a positive integer.');
    end

    [X, Y, N] = size(roiIn);

    % ---- Step 1: extract connected components across all input slices ----
    comps = struct('pixIdx', {}, 'vals', {}, 'centroid', {}, 'sources', {});
    idxComp = 0;

    for k = 1:N
        roiSlice = roiIn(:,:,k);
        bw = roiSlice > 0;
        cc = bwconncomp(bw);

        for j = 1:cc.NumObjects
            idxComp = idxComp + 1;
            pixIdx = cc.PixelIdxList{j};

            % values at those pixels (preserve original ROISlice values)
            vals = roiSlice(pixIdx);

            % centroid
            [r, c] = ind2sub([X, Y], pixIdx);
            centroid = [mean(r), mean(c)];

            comps(idxComp).pixIdx   = pixIdx(:);
            comps(idxComp).vals     = vals(:);
            comps(idxComp).centroid = centroid;
            comps(idxComp).sources  = k;  % start with this original slice index
        end
    end

    if isempty(comps)
        roiOut = zeros(X, Y, 0);
        mapping = cell(0,1);
        return;
    end

    % ---- Helper lambdas ----
    compSize = @(s) numel(s.pixIdx);

    function s = recomputeCentroid(s)
        [r, c] = ind2sub([X, Y], s.pixIdx);
        s.centroid = [mean(r), mean(c)];
    end

    function s = mergeTwo(a, b)
        % Merge components a and b into one.
        % If overlapping pixels exist, keep max value at overlap.
        % Most ROI pipelines avoid overlap, but we handle it safely.

        % Combine pixel indices
        allIdx = [a.pixIdx; b.pixIdx];
        allVal = [a.vals;   b.vals];

        % Resolve duplicates by max value
        % (sort and aggregate)
        [allIdxSorted, order] = sort(allIdx);
        allValSorted = allVal(order);

        % Find runs of equal indices
        dupe = [false; diff(allIdxSorted)==0];
        if any(dupe)
            % group by index, take max per group
            [uIdx, ~, g] = unique(allIdxSorted, 'stable');
            maxVal = accumarray(g, allValSorted, [], @max);
            pixIdx = uIdx;
            vals   = maxVal;
        else
            pixIdx = allIdxSorted;
            vals   = allValSorted;
        end

        s.pixIdx  = pixIdx;
        s.vals    = vals;
        s.sources = unique([a.sources(:); b.sources(:)]).';
        s = recomputeCentroid(s);
    end

    function nn = nearestNeighborIndex(i, compsLocal)
        % Find nearest neighbor of compsLocal(i) among all other components.
        ci = compsLocal(i).centroid;
        dBest = inf;
        nn = -1;
        for t = 1:numel(compsLocal)
            if t == i, continue; end
            ct = compsLocal(t).centroid;
            d = hypot(ci(1)-ct(1), ci(2)-ct(2));
            if d < dBest
                dBest = d;
                nn = t;
            end
        end
    end

    % ---- Step 2: iteratively merge small components ----
    % Keep merging the smallest ROI (< minSize) with its nearest neighbor.
    while true
        sizes = arrayfun(compSize, comps);
        smallIdx = find(sizes < minSize);

        if isempty(smallIdx)
            break; % all good
        end

        if numel(comps) == 1
            break; % cannot merge further
        end

        % pick the smallest component (or first among ties)
        [~, ord] = sort(sizes(smallIdx), 'ascend');
        i = smallIdx(ord(1));

        j = nearestNeighborIndex(i, comps);
        if j < 0
            break;
        end

        % Merge i and j; remove the higher index first to keep indexing valid
        newComp = mergeTwo(comps(i), comps(j));

        % Replace one slot with merged comp, delete the other
        keep = min(i, j);
        kill = max(i, j);

        comps(keep) = newComp;
        comps(kill) = [];
    end

    % ---- Step 3: pack output into X x Y x M ----
    M = numel(comps);
    roiOut = zeros(X, Y, M);
    mapping = cell(M, 1);

    for m = 1:M
        temp = zeros(X, Y);
        temp(comps(m).pixIdx) = comps(m).vals;
        roiOut(:,:,m) = temp;
        mapping{m} = comps(m).sources;
    end
end
