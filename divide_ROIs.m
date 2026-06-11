function roi_out = divide_ROIs(roi_stack, target_size, skeleton)
% divide_ROIs  Divide each ROI in a binary stack into equal-sized sub-ROIs.
%
% For each ROI slice, cuts the region into strips of ~target_size pixels
% oriented along the skeleton direction (or long axis if no skeleton).
% Strips smaller than target_size/2 are merged into their nearest neighbor.
%
% Usage:
%   roi_out = divide_ROIs(roi_stack, target_size)
%   roi_out = divide_ROIs(roi_stack, target_size, skeleton)
%
% Inputs:
%   roi_stack   - [H x W x N] logical array. Each slice is one ROI binary mask.
%   target_size - Scalar. Target number of pixels per output sub-ROI.
%   skeleton    - (optional) [H x W] binary skeleton image. When provided,
%                 the cutting direction per ROI is taken from the local
%                 skeleton orientation within that ROI. If absent or if the
%                 skeleton has no pixels in a given ROI, the long axis of
%                 the ROI (via PCA on pixel coords) is used instead.
%
% Output:
%   roi_out     - [H x W x M] logical array. M >= N. Each slice is one
%                 sub-ROI. Sub-ROIs from the same input ROI do not overlap.
%                 Sub-ROIs from different input ROIs also do not overlap.
%
% Algorithm (per input ROI):
%   1. Find the cutting axis: fit a line to skeleton pixels inside the ROI
%      (or use PCA on all ROI pixels if no skeleton).
%   2. Project all ROI pixels onto that axis to get a 1D position per pixel.
%   3. Sort pixels by projection, divide into strips of target_size, assign
%      each pixel a strip index.
%   4. Merge any strip smaller than target_size/2 into the adjacent strip
%      with fewer pixels (to keep sizes balanced).
%   5. Output one mask per final strip.
%
% Example:
%   roi_out = divide_ROIs(roi_stack, 200, skeleton);
%   label_img = sum(roi_out .* permute(1:size(roi_out,3),[1 3 2]), 3);
%   figure; imagesc(label_img); axis image; colormap(lines);
%
% 2025.05  Byung Hun Lee

%% ---- 0. Input handling ----
roi_stack   = logical(roi_stack);
[img_h, img_w, nROIs] = size(roi_stack);
target_size = double(target_size);
merge_threshold = target_size / 2;   % strips smaller than this get merged

use_skeleton = nargin >= 3 && ~isempty(skeleton) && any(skeleton(:));
if use_skeleton
    skeleton = logical(skeleton);
end

all_masks = {};

%% ---- 1. Process each input ROI ----
for iROI = 1:nROIs

    % Pixel coordinates of this ROI (col=x, row=y)
    [px_rows, px_cols] = find(roi_stack(:,:,iROI));
    nPix = numel(px_rows);

    if nPix == 0
        continue
    end

    if nPix <= target_size
        % ROI already smaller than target — keep as one piece
        all_masks{end+1} = roi_stack(:,:,iROI); %#ok<AGROW>
        continue
    end

    %% ---- 2. Find cutting axis ----
    % Try skeleton pixels inside this ROI first
    axis_vec = [];
    if use_skeleton
        skel_in_roi = skeleton & roi_stack(:,:,iROI);
        [sk_rows, sk_cols] = find(skel_in_roi);
        if numel(sk_rows) >= 2
            % PCA on skeleton pixels to get primary orientation
            sk_pts   = [sk_cols, sk_rows] - mean([sk_cols, sk_rows], 1);
            [~, ~, V] = svd(sk_pts, 'econ');
            axis_vec = V(:,1)';   % [1 x 2] unit vector along skeleton
        end
    end

    % Fall back to PCA on all ROI pixels
    if isempty(axis_vec)
        roi_pts  = [px_cols, px_rows] - mean([px_cols, px_rows], 1);
        [~, ~, V] = svd(roi_pts, 'econ');
        axis_vec = V(:,1)';
    end

    %% ---- 3. Project ROI pixels onto cutting axis, sort, assign strips ----
    px_xy       = [px_cols, px_rows];                          % [nPix x 2]
    px_xy_c     = px_xy - mean(px_xy, 1);                      % centre
    projection  = px_xy_c * axis_vec';                         % [nPix x 1] scalar position
    [~, sort_idx] = sort(projection);                          % sort along axis

    n_strips  = max(1, round(nPix / target_size));
    strip_idx = zeros(nPix, 1);
    boundaries = round(linspace(0, nPix, n_strips + 1));
    for iStrip = 1:n_strips
        pix_in_strip = sort_idx(boundaries(iStrip)+1 : boundaries(iStrip+1));
        strip_idx(pix_in_strip) = iStrip;
    end

    %% ---- 4. Merge strips smaller than merge_threshold ----
    strip_idx = merge_small_strips(strip_idx, n_strips, merge_threshold);
    n_strips_final = max(strip_idx);

    %% ---- 5. Build one mask per final strip ----
    for iStrip = 1:n_strips_final
        in_strip = strip_idx == iStrip;
        if ~any(in_strip)
            continue
        end
        mask = false(img_h, img_w);
        lin  = sub2ind([img_h, img_w], px_rows(in_strip), px_cols(in_strip));
        mask(lin) = true;
        all_masks{end+1} = mask; %#ok<AGROW>
    end
end

%% ---- 6. Stack output ----
nOut = numel(all_masks);
if nOut == 0
    warning('divide_ROIs: no output ROIs generated.');
    roi_out = false(img_h, img_w, 0);
    return
end

roi_out = false(img_h, img_w, nOut);
for n = 1:nOut
    roi_out(:,:,n) = all_masks{n};
end

fprintf('divide_ROIs: %d input ROIs -> %d output sub-ROIs (target_size = %d px)\n', ...
    nROIs, nOut, target_size);
end


%% ---- Helper: merge strips smaller than threshold into nearest neighbour ----
function strip_idx = merge_small_strips(strip_idx, n_strips, threshold)
% Iteratively find the smallest strip; if it is below threshold, absorb it
% into the adjacent strip (left or right) that currently has fewer pixels.
% Repeat until all strips meet the threshold or only one strip remains.

while true
    strip_sizes = arrayfun(@(s) sum(strip_idx == s), 1:n_strips);
    [min_size, min_s] = min(strip_sizes);

    if min_size >= threshold || n_strips == 1
        break
    end

    % Determine neighbours
    left  = min_s - 1;
    right = min_s + 1;

    if left < 1
        target_strip = right;
    elseif right > n_strips
        target_strip = left;
    else
        % Merge into whichever neighbour is smaller (keep sizes balanced)
        if strip_sizes(left) <= strip_sizes(right)
            target_strip = left;
        else
            target_strip = right;
        end
    end

    % Absorb min_s into target_strip
    strip_idx(strip_idx == min_s) = target_strip;

    % Renumber strips to be contiguous 1..n_strips-1
    old_labels = unique(strip_idx);
    new_labels  = 1:numel(old_labels);
    new_idx     = zeros(size(strip_idx));
    for k = 1:numel(old_labels)
        new_idx(strip_idx == old_labels(k)) = new_labels(k);
    end
    strip_idx = new_idx;
    n_strips  = n_strips - 1;
end
end
