function [ match_data, match_img ] = match_pattern( img, patt, CORR_THRESH, MIN_DIST, FFTimg, FFTpatt)
% Detects a single pattern 'patt' withing the image 'img'. Only matches above the
% correlation threshold threshold 'CORR_THRESH' are taken as valid. Each
% returned match is the best candidate within a box window with edge
% length 2*MIN_DIST+1 centered around the candidate.
%
% USAGE: [ match_data, match_img ] = match_pattern(img, patt, CORR_THRESH, MIN_DIST)
%
% Input:
%    img         - The image to search in
%    patt        - The pattern to search for
%    CORR_THRESH - Minimum (absolute) normalized correlation value for a
%                  match to be valid. Valid values are between 0 and 1.
%    MIN_DIST    - Only the best candidates within a box window with edge
%                  length 2*MIN_DIST+1 centered around the candidate are taken.
%    FFTimg      - (optional) fft2 of img if known. see normxcorr2_preFFT for more info.
%    FFTpatt     - (optional) fft2 of pattern if known. see normxcorr2_preFFT for more info.
%
% Output:
%    match_data  - Mx3 matrix, where M is the number of valid matches. The
%                  columns correspond to (rowShift, colShift, correlation 
%                  value). For rowShift = colShift = 0 the upper left pixel
%                  of the pattern lies ontop the upper left pixel of the
%                  image. The pattern centers can be calculated as
%                  (rowShift + pattRows/2, colShift + pattCols/2).
%    match_img   - Image of size (imgRows-pattRows)x(imgCols-pattCols) 
%                  showing pattern shifts with pixels intensities 
%                  corresponding to the correlation values.
%    
% Note: Only patterns fully embedded in the image are detected.

%Note: CORR_THRESH depends on 'empty' space around the pattern.
% For larger zero padding and the presence of noise in the image
% this should be reduced, as the normalized correlation drops down.
% That said, zero padding seems to improve the performance of the detection
% (while sacrificing the usable border).

% Author: Simon Christoph Stein
% Date: 2014
% EMail: scstein@phys.uni-goettingen.de

MIN_DIST = max(1, MIN_DIST); % Distance should be at least 1.


imgRows = size(img,1);
imgCols = size(img,2);

pattRows = size(patt,1);
pattCols = size(patt,2);

% -- Pattern detection by cross-correlation --
% NOTE: Normalized cross correlation should to be used for feature matching,
% as otherwise the amplitude of the measured image will influence the
% result. For example: a spot much brighter than the pattern in the image
% will always lead to a local maximum.

% CC = xcorr2(img-mean(img(:)), patt);
% CC = normxcorr2(patt,img); % normxcorr2 requires a smaller template than the image!
if nargin < 5; FFTimg = []; end;
if nargin < 6; FFTpatt = []; end;
CC = normxcorr2_preFFT(patt,img, FFTimg, FFTpatt); % normxcorr2 with precomputed ffts
CC = CC(pattRows:imgRows, pattCols:imgCols); % cut out valid part of the cross-correlation (pattern fully embedded)

% Set all points below the thresold to the minimum
mask_overTOL = abs(CC)>CORR_THRESH;


% -- Detect local maxima within sliding window --
localmax_idx = localMax(CC, MIN_DIST);
mask_localmax = zeros(size(CC));
mask_localmax(localmax_idx) = 1;

mask = mask_localmax & mask_overTOL;


% Gather the positions of all local maxima
max_ind = find(mask);
[max_rPos, max_cPos] = ind2sub(size(CC), max_ind);
max_shifts = [max_rPos, max_cPos] -1; % (rowShift, colShift) of maxima

match_data = [max_shifts, CC(mask)]; % (rowShift, colShift, corrVal)


% image of matches with their correlation value
match_img = CC;
match_img(~mask) = 0;
end



function max_indices = localMax(img, MIN_DIST)
% Returns indices of highest local maximum within a box window of edge
% length 2*MIN_DIST+1 around each pixel. Connected pixels with equal values
% are considered maxima as well.
    neigh = ones(2*MIN_DIST+1); % Neighborhood matrix
    img_dilated = imdilate(img, neigh);
    max_indices = find(img == img_dilated);
end

