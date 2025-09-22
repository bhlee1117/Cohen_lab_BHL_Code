function [full_out, scaleImgs_all] = SeeResiduals_tiling(mov, refSig, remOffset, segmentLength, overlapRatio)
% SeeResiduals_tiling - Applies SeeResiduals in overlapping temporal tiles.
%
% Inputs:
%   mov           : [Y, X, T] movie
%   refSig        : [N, T] reference traces
%   remOffset     : 1 or 0, whether to remove DC offset
%   segmentLength : number of frames in each tile
%   overlapRatio  : fraction (0-1) of overlap between tiles
%
% Outputs:
%   full_out      : residual movie after stitched tiles
%   scaleImgs_all : scale images for each segment (cell array)

if nargin < 5
    overlapRatio = 0.05;
end
if nargin < 4
    segmentLength = 500;
end
if nargin < 3
    remOffset = 1;
end

if size(refSig,1)>size(refSig,2)
    refSig=refSig';
end

[Y, X, T] = size(mov);
step = round(segmentLength * (1 - overlapRatio));
tileStart = 1:step:T;

full_out = zeros(Y, X, T);
weight = zeros(1, 1, T);  % keep track of contribution per frame
scaleImgs_all = {};

for i = 1:length(tileStart)
    t0 = tileStart(i);
    t1 = min(t0 + segmentLength - 1, T);
    idx = t0:t1;

    mov_seg = mov(:, :, idx);
    ref_seg = refSig(:, idx);

    % pad if last segment is shorter
    if size(ref_seg, 2) < segmentLength
        padN = segmentLength - size(ref_seg, 2);
        ref_seg(:, end+1:end+padN) = repmat(ref_seg(:, end), 1, padN);
        mov_seg(:, :, end+1:end+padN) = repmat(mov_seg(:, :, end), 1, 1, padN);
    end

    [resid, scaleImgs] = SeeResiduals(mov_seg, ref_seg, remOffset);
    resid = resid(:, :, 1:length(idx));  % remove padding
    full_out(:, :, idx) = full_out(:, :, idx) + resid;
    weight(:, :, idx) = weight(:, :, idx) + 1;

    scaleImgs_all{end+1} = scaleImgs;
end

% Normalize by overlap
weight(weight == 0) = 1;
full_out = full_out ./ weight;

end
