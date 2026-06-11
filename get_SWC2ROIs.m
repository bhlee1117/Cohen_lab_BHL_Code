function roi_stack = get_SWC2ROIs(skeleton, outline, roi_size)
% get_SWC2ROIs  Divide a neuron outline into ROIs along a skeleton.
%
% Traces each branch of a binary skeleton, divides it into arc-length windows
% of roi_size, then assigns outline pixels to the nearest skeleton node within
% each window. Because branches are traced independently, ROIs never bleed
% across unconnected branches even when they are spatially close.
%
% Usage:
%   roi_stack = get_SWC2ROIs(skeleton, outline, roi_size)
%
% Inputs:
%   skeleton  - [H x W] binary image (logical or 0/1). The thinned skeleton
%               of the neuron. Should be 1-pixel wide (e.g. output of bwmorph
%               with 'thin' or 'skel'). If not already thinned, this function
%               thins it automatically.
%   outline   - [H x W] binary image (logical or 0/1). The filled or outline
%               mask of the neuron whose pixels will be partitioned into ROIs.
%   roi_size  - Scalar. Target arc-length of each ROI along the skeleton
%               (pixels). Controls how many ROIs are placed per branch.
%
% Output:
%   roi_stack - [H x W x N] logical array. Each slice (:,:,n) is the binary
%               mask of one ROI. N = total ROIs across all branches.
%               Slices do not overlap. Outline pixels farther from the
%               skeleton than roi_size may be left unassigned.
%
% Algorithm:
%   1. Thin skeleton to 1px width, find branch points and endpoints.
%   2. Split skeleton into individual branch segments by removing branch pts.
%   3. For each segment: order pixels along the path, compute arc length,
%      divide into roi_size windows.
%   4. Assign each outline pixel to its nearest skeleton pixel globally
%      (Voronoi ownership), then group by window membership.
%   5. Stack all non-empty ROI masks into [H x W x N].
%
% Example:
%   skel    = bwmorph(neuron_mask, 'thin', Inf);
%   outline = neuron_mask;
%   roi_stack = get_SWC2ROIs(skel, outline, 20);
%   label_img = sum(roi_stack .* permute(1:size(roi_stack,3),[1 3 2]), 3);
%   figure; imagesc(label_img); axis image; colorbar;
%
% 2025.05  Byung Hun Lee

%% ---- 0. Input handling ----
skeleton = logical(skeleton);
outline  = logical(outline);
[img_h, img_w] = size(skeleton);

if ~isequal(size(skeleton), size(outline))
    error('get_SWC2ROIs: skeleton and outline must be the same size.');
end

%% ---- 1. Ensure skeleton is 1px wide ----
skeleton = bwmorph(skeleton, 'thin', Inf);

%% ---- 2. Find branch points and endpoints ----
% Branch points: skeleton pixels with 3+ neighbours
% Endpoints:     skeleton pixels with exactly 1 neighbour
branch_pts = bwmorph(skeleton, 'branchpoints');
end_pts    = bwmorph(skeleton, 'endpoints');

%% ---- 3. Split skeleton into individual branch segments ----
% Remove branch points to disconnect branches, then label connected pieces
skeleton_split = skeleton & ~branch_pts;
[seg_label, nSegments] = bwlabel(skeleton_split, 8);

%% ---- 4. Pre-compute Voronoi ownership: nearest skeleton pixel per outline pixel ----
% All skeleton pixel coordinates
[skel_rows, skel_cols] = find(skeleton);
skel_xy = [skel_cols, skel_rows];   % [nSkelPix x 2] (x=col, y=row)
nSkelPix = size(skel_xy, 1);

% All outline pixel coordinates
[out_rows, out_cols] = find(outline);
out_xy = [out_cols, out_rows];      % [nOutPix x 2]
nOutPix = size(out_xy, 1);

if nOutPix == 0
    warning('get_SWC2ROIs: outline is empty. Returning empty stack.');
    roi_stack = false(img_h, img_w, 0);
    return
end

% For each outline pixel find its nearest skeleton pixel (global Voronoi)
if exist('knnsearch', 'file') == 2
    nearest_skel_idx = knnsearch(skel_xy, out_xy, 'K', 1);
else
    D = pdist2(out_xy, skel_xy);
    [~, nearest_skel_idx] = min(D, [], 2);
end

% Convert nearest skeleton pixel index to its segment label
nearest_seg_label = seg_label(sub2ind([img_h, img_w], skel_rows(nearest_skel_idx), skel_cols(nearest_skel_idx)));

%% ---- 5. For each segment: order pixels, divide into windows, build ROI masks ----
all_roi_masks = {};

for iSeg = 1:nSegments

    % Pixel coordinates of this segment
    [s_rows, s_cols] = find(seg_label == iSeg);
    seg_xy = [s_cols, s_rows];   % [k x 2] (x, y)
    k = size(seg_xy, 1);

    if k < 2
        continue
    end

    % -- Order segment pixels along the path --
    % Start from whichever pixel in this segment is closest to
    % any endpoint or branch point (i.e. one tip of the segment)
    seg_lin_idx = sub2ind([img_h, img_w], s_rows, s_cols);

    [ep_rows, ep_cols] = find(end_pts | branch_pts);
    ep_xy = [ep_cols, ep_rows];   % (x=col, y=row)

    if ~isempty(ep_xy)
        D_tip = pdist2(seg_xy, ep_xy);
        [~, start_local] = min(min(D_tip, [], 2));
    else
        start_local = 1;
    end

    ordered = zeros(k, 2);
    visited = false(k, 1);
    ordered(1, :) = seg_xy(start_local, :);
    visited(start_local) = true;

    for iStep = 2:k
        prev = ordered(iStep-1, :);
        d = sqrt(sum((seg_xy - prev).^2, 2));
        d(visited) = Inf;
        [~, next_local] = min(d);
        ordered(iStep, :) = seg_xy(next_local, :);
        visited(next_local) = true;
    end

    % -- Cumulative arc length along ordered segment --
    step_d = sqrt(sum(diff(ordered, 1, 1).^2, 2));   % [k-1 x 1]
    arc_s  = [0; cumsum(step_d)];                     % [k x 1]
    total_s = arc_s(end);

    % -- Divide into windows --
    n_windows = max(1, floor(total_s / roi_size));
    win_edges = linspace(0, total_s, n_windows + 1);

    % Linear indices of ordered skeleton pixels (for membership lookup)
    ordered_lin = sub2ind([img_h, img_w], ordered(:,2), ordered(:,1));

    for iWin = 1:n_windows
        s_lo = win_edges(iWin);
        s_hi = win_edges(iWin + 1);

        in_window     = (arc_s >= s_lo) & (arc_s <= s_hi);
        win_skel_lin  = ordered_lin(in_window);   % skeleton pixel lin-indices in this window

        % Convert win_skel_lin to indices into the global skel_rows/skel_cols arrays
        % so we can match against nearest_skel_idx
        [~, win_skel_global_idx] = ismember(win_skel_lin, seg_lin_idx(1:end));
        % Re-map: find positions in global skel array
        win_global = find(ismember(sub2ind([img_h img_w], skel_rows, skel_cols), win_skel_lin));

        % Outline pixels whose nearest skeleton pixel falls in this window
        in_roi = ismember(nearest_skel_idx, win_global) & (nearest_seg_label == iSeg);

        if ~any(in_roi)
            continue
        end

        % Build binary mask
        roi_mask = false(img_h, img_w);
        lin_idx  = sub2ind([img_h, img_w], out_rows(in_roi), out_cols(in_roi));
        roi_mask(lin_idx) = true;
        all_roi_masks{end+1} = roi_mask; %#ok<AGROW>
    end
end

%% ---- 6. Stack into [H x W x N] output ----
nROIs = numel(all_roi_masks);
if nROIs == 0
    warning('get_SWC2ROIs: no ROIs generated. Check that skeleton and outline overlap.');
    roi_stack = false(img_h, img_w, 0);
    return
end

roi_stack = false(img_h, img_w, nROIs);
for n = 1:nROIs
    roi_stack(:, :, n) = all_roi_masks{n};
end

fprintf('get_SWC2ROIs: %d segments -> %d ROIs (roi_size = %.1f px)\n', nSegments, nROIs, roi_size);
end
