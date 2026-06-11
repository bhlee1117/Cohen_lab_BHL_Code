function psd_mov = mov2psdmov(mov, psd, frequencies, freq_bins)
% mov2psdmov  Convert a movie and its per-pixel PSD into a power spectrum movie.
%
% For each pixel, integrates (sums) PSD power within each frequency bin,
% producing a movie where each frame corresponds to one frequency band.
%
% Usage:
%   psd_mov = mov2psdmov(mov, psd, frequencies, freq_bins)
%
% Inputs:
%   mov         - [H x W x T] numeric array. The original movie (used only
%                 to recover spatial dimensions H and W).
%   psd         - [H*W x nFreq] numeric array. Per-pixel PSD as returned by
%                 nanPSD (first output). Rows are pixels in column-major order
%                 matching reshape(mov, [], T). Columns are frequency bins
%                 corresponding to the frequencies vector.
%   frequencies - [1 x nFreq] frequency vector in Hz, as returned by nanPSD
%                 (second output).
%   freq_bins   - Frequency bin definition. Two formats accepted:
%                   [nBins x 2]  : each row is [f_low, f_high] (Hz). Power
%                                  is integrated over f_low <= f <= f_high.
%                   [1 x (nB+1)]: bin edges vector. Bins are defined as
%                                  consecutive pairs of edges (like histcounts).
%                                  Produces nB = numel(freq_bins)-1 output frames.
%
% Output:
%   psd_mov     - [H x W x nBins] single array. Each slice (:,:,b) is the
%                 spatially-resolved integrated power in frequency bin b.
%                 Units match those of the input psd (e.g. V²/Hz integrated
%                 over Hz = V²). NaN pixels in psd propagate as NaN.
%
% Example:
%   [psd, freq] = nanPSD(reshape(mov,[],T)', fs, 3000);
%   psd = psd';   % nanPSD returns [nRows x nFreq]; rows=pixels here
%
%   % Option 1: explicit band edges [nBins x 2]
%   freq_bins = [1 4;  4 8;  8 12;  12 30;  30 80];   % delta,theta,alpha,beta,gamma
%   psd_mov = mov2psdmov(mov, psd, freq, freq_bins);
%
%   % Option 2: edge vector
%   freq_bins = [1 4 8 12 30 80];
%   psd_mov = mov2psdmov(mov, psd, freq, freq_bins);
%
%   % Visualise band 3 (alpha)
%   imagesc(psd_mov(:,:,3)); colorbar; axis image;
%
% See also: nanPSD
%
% 2025.05  Byung Hun Lee

%% ---- 0. Input handling ----
[img_h, img_w, ~] = size(mov);
nPix              = img_h * img_w;
frequencies       = frequencies(:)';        % ensure row vector [1 x nFreq]
psd               = double(psd);

if size(psd, 1) ~= nPix
    error('mov2psdmov: psd must have H*W = %d rows, but has %d.', nPix, size(psd,1));
end
if size(psd, 2) ~= numel(frequencies)
    error('mov2psdmov: psd columns (%d) must match numel(frequencies) (%d).', ...
        size(psd,2), numel(frequencies));
end

%% ---- 1. Parse freq_bins into [nBins x 2] low/high pairs ----
if isvector(freq_bins) && numel(freq_bins) >= 2
    % Edge vector -> consecutive pairs
    freq_bins = freq_bins(:);
    freq_bins = [freq_bins(1:end-1), freq_bins(2:end)];   % [nBins x 2]
elseif size(freq_bins, 2) == 2
    % Already [nBins x 2]
    freq_bins = double(freq_bins);
else
    error('mov2psdmov: freq_bins must be [nBins x 2] or a [1 x (nBins+1)] edge vector.');
end
nBins = size(freq_bins, 1);

%% ---- 2. Pre-compute frequency mask for each bin ----
% freq_mask{b} is a logical [1 x nFreq] selecting frequencies in bin b
freq_masks = false(nBins, numel(frequencies));
for iBin = 1:nBins
    freq_masks(iBin, :) = (frequencies >= freq_bins(iBin,1)) & ...
                          (frequencies <= freq_bins(iBin,2));
    if ~any(freq_masks(iBin,:))
        warning('mov2psdmov: bin %d [%.2f %.2f Hz] contains no frequency samples.', ...
            iBin, freq_bins(iBin,1), freq_bins(iBin,2));
    end
end

%% ---- 3. Integrate PSD within each bin for all pixels at once ----
% psd is [nPix x nFreq]; freq_masks is [nBins x nFreq]
% bin_power = psd * freq_masks' -> [nPix x nBins]
% (equivalent to summing psd columns in each bin)
bin_power = psd * freq_masks';    % [nPix x nBins]

%% ---- 4. Reshape to [H x W x nBins] ----
psd_mov = single(reshape(bin_power, img_h, img_w, nBins));

fprintf('mov2psdmov: %d bins computed over %d x %d pixels\n', nBins, img_h, img_w);
end
