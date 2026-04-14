function [dist, skeleton_path] = geodesic_distance(bwimg, point1, point2)
% geodesic_distance  Geodesic distance and path between two points on a
%                    binary skeleton/mask.
%
% Usage:
%   [dist, skeleton_path] = geodesic_distance(bwimg, point1, point2)
%
% Inputs:
%   bwimg  - binary image (logical or 0/1), the constrained path domain
%   point1 - [x y] centroid of ROI 1  (col, row)
%   point2 - [x y] centroid of ROI 2  (col, row)
%
% Outputs:
%   dist          - scalar geodesic distance (pixels) along bwimg
%   skeleton_path - binary image of the shortest path
%
% Strategy (no while-loops, guaranteed termination):
%   1. Snap each point to the nearest ON pixel of bwimg via bwdist.
%   2. Draw a thin straight-line bridge from each snapped point to its
%      original centroid (handles points that lie off the skeleton).
%   3. If the two snapped points still belong to different connected
%      components, draw one more straight-line bridge between them through
%      the mask (this always terminates in O(1) operations).
%   4. Run bwdistgeodesic on the unified mask.

    point1 = round(point1(:))';   % ensure 1x2  [x y]
    point2 = round(point2(:))';

    [H, W] = size(bwimg);

    %-- Clamp to image bounds
    point1 = max(min(point1, [W H]), [1 1]);
    point2 = max(min(point2, [W H]), [1 1]);

    bwimg  = logical(bwimg);

    %-- Step 1: snap each point to the nearest ON pixel -----------------
    [D_all, idx_all] = bwdist(bwimg);   % distance from every pixel to nearest ON pixel

    p1_snapped = snap_to_mask(point1, bwimg, D_all, idx_all, H, W);
    p2_snapped = snap_to_mask(point2, bwimg, D_all, idx_all, H, W);

    %-- Step 2: draw straight-line bridges from centroid -> snapped point
    bwimg = draw_line(bwimg, point1,     p1_snapped, H, W);
    bwimg = draw_line(bwimg, point2,     p2_snapped, H, W);

    %-- Step 3: connect the two snapped points if in different components
    bwss = bwlabel(bwimg);
    if bwss(p1_snapped(2), p1_snapped(1)) ~= bwss(p2_snapped(2), p2_snapped(1))
        % Draw one straight bridge between the two snapped points.
        % This always terminates (single O(n) line draw, no loop).
        bwimg = draw_line(bwimg, p1_snapped, p2_snapped, H, W);
    end

    %-- Step 4: geodesic distance on unified mask -----------------------
    bwimg = logical(bwimg);
    D1 = bwdistgeodesic(bwimg, point1(1), point1(2), 'quasi-euclidean');
    D2 = bwdistgeodesic(bwimg, point2(1), point2(2), 'quasi-euclidean');

    D = D1 + D2;
    D = round(D * 10) / 10;
    D(isnan(D)) = inf;

    %-- Extract skeleton path between the two points --------------------
    skeleton_path = imregionalmin(D);
    skeleton_path = bwmorph(skeleton_path, 'fill',  Inf);
    skeleton_path = bwmorph(skeleton_path, 'close', Inf);
    skeleton_path = bwmorph(skeleton_path, 'skel',  Inf);

    dist = median(D(skeleton_path), 'omitnan');
end

%% =========================================================================
function p_snapped = snap_to_mask(p, bwimg, ~, idx_all, H, W)
% Return the nearest ON pixel to p. If p is already ON, return p unchanged.
    col = p(1);  row = p(2);
    if bwimg(row, col)
        p_snapped = p;
    else
        lin = idx_all(row, col);          % linear index of nearest ON pixel
        [r2, c2] = ind2sub([H W], lin);
        p_snapped = [c2, r2];             % [x y]
    end
end

%% =========================================================================
function bwimg = draw_line(bwimg, p1, p2, H, W)
% Rasterise a straight line between p1 and p2 (Bresenham) and set those
% pixels ON in bwimg. Guaranteed O(max(dx,dy)) — no loop termination risk.
    x0 = p1(1);  y0 = p1(2);
    x1 = p2(1);  y1 = p2(2);

    dx = abs(x1 - x0);   dy = abs(y1 - y0);
    sx = sign(x1 - x0);  sy = sign(y1 - y0);

    % Pre-allocate pixel list (upper bound: max(dx,dy)+1 pixels)
    n   = max(dx, dy) + 1;
    xs  = zeros(n, 1, 'int32');
    ys  = zeros(n, 1, 'int32');
    err = dx - dy;
    xi  = x0;  yi  = y0;

    for k = 1:n
        xs(k) = xi;
        ys(k) = yi;
        if xi == x1 && yi == y1, break; end
        e2 = 2 * err;
        if e2 > -dy
            err = err - dy;
            xi  = xi  + sx;
        end
        if e2 <  dx
            err = err + dx;
            yi  = yi  + sy;
        end
    end

    % Clamp and write
    xs = min(max(xs(1:k), 1), W);
    ys = min(max(ys(1:k), 1), H);
    idx = sub2ind([H W], double(ys), double(xs));
    bwimg(idx) = true;
end
