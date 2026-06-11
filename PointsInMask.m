function [inTrue, xy_inTrue] = PointsInMask(xy, bw)
% POINTSINMASK  Find which XY coordinates fall on 'true' pixels in a binary mask.
%
% Usage:
%   [inTrue, xy_inTrue] = PointsInMask(xy, bw)
%
% Inputs:
%   xy  - N x 2 matrix of coordinates [x, y] (x = column, y = row)
%         Points outside the image boundary are treated as false.
%   bw  - 2D binary image (logical or numeric)
%
% Outputs:
%   inTrue    - N x 1 logical vector; true where the point lands on a true pixel
%   xy_inTrue - M x 2 matrix of coordinates that land on true pixels

sz = size(bw);

% Filter points within image bounds
valid = xy(:,1) >= 1 & xy(:,1) <= sz(2) & ...
        xy(:,2) >= 1 & xy(:,2) <= sz(1);

% Convert to linear indices and check mask value
idx = sub2ind(sz, round(xy(valid,2)), round(xy(valid,1)));

inTrue = false(size(xy,1), 1);
inTrue(valid) = logical(bw(idx));

xy_inTrue = xy(inTrue, :);

end
