function [tr_sub, tr_sub2] = get_subthreshold(tr, sp, dilate, avgwnd)
% tr    : N x T voltage trace
% sp    : N x T binary spike trace
% dilate: peri-spike frames to omit (e.g. 7 → ±3 frames)
% avgwnd: smoothing window after interpolation

se = strel('square', dilate);

binary_dilation = double(any(sp, 1));
binary_dilation = imdilate(binary_dilation, se);

valid_point = find(~binary_dilation);

tr_sub = NaN(size(tr));
if numel(valid_point) >= 2
    for n = 1:size(tr, 1)
        tr_sub(n, :) = interp1(valid_point, tr(n, valid_point)', 1:size(tr, 2), 'linear');
    end
elseif numel(valid_point) == 1
    % Only one non-masked frame: interp1 needs >=2 points, so fall back to a
    % constant baseline (that frame's value).
    tr_sub = repmat(tr(:, valid_point), 1, size(tr, 2));
else
    % Every frame is masked by the spike/blue dilation (nothing to interpolate
    % from). Leave tr_sub as NaN so the caller can detect/skip this trace.
    warning('get_subthreshold:noValidPoints', ...
        'All %d frames masked by spike/blue dilation; returning NaN subthreshold.', size(tr, 2));
end

tr_sub2 = tr_sub;
tr_sub  = movmean(tr_sub, avgwnd, 2, 'omitnan');

end