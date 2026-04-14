function Bout = imgaussfiltnan(A, sigma, varargin)
%IMGAUSSFILTNAN Gaussian filter that ignores NaN values.
%
%   B = IMGAUSSFILTNAN(A, sigma) applies a Gaussian filter with standard
%   deviation sigma to array A while ignoring NaNs. Works like IMGAUSSFILT
%   but treats NaNs as missing data rather than zeros.
%
%   Additional name/value pairs accepted by IMGAUSSFILT can be passed after
%   sigma (e.g. 'FilterSize',[5 5]).
%
%   When NaN pixels are present the function crops to the tight bounding box
%   of valid pixels (plus a half-kernel periphery), filters only that
%   sub-region, and pastes the result back — avoiding redundant computation
%   on empty pixels.
%
% Example:
%   A = peaks(100); A(30:40,30:40) = NaN;
%   B = imgaussfiltnan(A, 2);
%
% Byung-Hun / Pawgoomon, 2025

%-- Build Gaussian kernel
if isempty(varargin)
    ksize = 2*ceil(2*sigma) + 1;
else
    p = inputParser;
    addOptional(p, 'FilterSize', 2*ceil(2*sigma)+1);
    parse(p, varargin{:});
    ksize = p.Results.FilterSize;
end
G = fspecial('gaussian', ksize, sigma);

[H, W, nZ] = size(A);
Bout = NaN(size(A));

%-- Bounding box of non-NaN pixels (shared across all z-slices for a fixed mask)
valid_pix = any(~isnan(A), 3);      % H x W: true wherever any slice has data
has_nan   = ~all(valid_pix(:));

if has_nan
    row_any = any(valid_pix, 2);    % H x 1
    col_any = any(valid_pix, 1);    % 1 x W
    pad     = ceil(max(ksize) / 2); % periphery = half kernel so edges are fully covered
    r1 = max(find(row_any, 1, 'first') - pad, 1);
    r2 = min(find(row_any, 1, 'last')  + pad, H);
    c1 = max(find(col_any, 1, 'first') - pad, 1);
    c2 = min(find(col_any, 1, 'last')  + pad, W);
else
    r1 = 1;  r2 = H;
    c1 = 1;  c2 = W;
end

%-- Filter only the cropped region and paste back
A_crop    = A(r1:r2, c1:c2, :);
Bout_crop = NaN(size(A_crop));

for z = 1:nZ
    M       = ~isnan(A_crop(:,:,z));
    Afilled = A_crop(:,:,z);
    Afilled(~M) = 0;                % replace NaNs with zero for convolution

    num = imfilter(Afilled,    G, 'same', 'replicate');
    den = imfilter(double(M),  G, 'same', 'replicate');

    B          = num ./ den;
    B(den==0)  = NaN;               % fully-NaN windows stay NaN
    Bout_crop(:,:,z) = B;
end

Bout(r1:r2, c1:c2, :) = Bout_crop;
end
