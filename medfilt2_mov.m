function mov_filt = medfilt2_mov(mov, sz)
% medfilt2_mov  Apply 2D median filter (NaN-aware) to each frame of a movie.
%
% If the movie contains NaN pixels (e.g. a masked ROI), the function
% automatically crops to the tight bounding box of valid pixels (plus a
% half-kernel periphery), filters only that sub-region, and pastes the
% result back — avoiding redundant computation on empty pixels.
%
% Progress messages are shown only when the crop area exceeds 10000 pixels.

if nargin < 2
    sz = [3 3];
end

mov_filt  = zeros(size(mov));
[H, W, T] = size(mov);

%-- Find bounding box of non-NaN pixels (same across all frames for a fixed mask)
valid_pix = any(~isnan(mov), 3);   % H x W logical: any frame has data here
has_nan   = ~all(valid_pix(:));    % true if there are NaN pixels at all

if has_nan
    row_any = any(valid_pix, 2);   % H x 1
    col_any = any(valid_pix, 1);   % 1 x W
    pad     = ceil(max(sz) / 2);   % periphery = half kernel
    r1 = max(find(row_any, 1, 'first') - pad, 1);
    r2 = min(find(row_any, 1, 'last')  + pad, H);
    c1 = max(find(col_any, 1, 'first') - pad, 1);
    c2 = min(find(col_any, 1, 'last')  + pad, W);
else
    r1 = 1;  r2 = H;
    c1 = 1;  c2 = W;
end

%-- Filter only the cropped region
mov_crop  = mov(r1:r2, c1:c2, :);
crop_filt = zeros(size(mov_crop));
crop_area = (r2-r1+1) * (c2-c1+1);
verbose   = crop_area > 10000;

for i = 1:T
    if verbose && mod(i, max(round(T/5), 1)) == 0
        fprintf('Medfilt processing... %5.1f%%\r', i/T*100);
    end
    crop_filt(:,:,i) = medfilt2nan(mov_crop(:,:,i), sz);
end

if verbose
    disp('Medfilt done.');
end

%-- Paste back into full-size output
mov_filt(r1:r2, c1:c2, :) = crop_filt;
end
