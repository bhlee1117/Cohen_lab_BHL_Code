function [psd, frequencies, psd_norm, psd_db] = nanPSD(x, fs, windowsize)
% NANPSD Compute PSD handling NaNs using averaged Lomb-Scargle on moving windows.
% [PSD, FREQUENCIES, PSD_NORM, PSD_DB] = NANPSD(X, FS)
% Computes PSD for each row in X, skipping all-NaN rows.
% Averages PSDs from overlapping segments (3000 frames window, 1500 overlap) for reduced noise.
% Uses plomb for uneven sampling due to NaNs removal.
% Normalizes and computes dB version.
% Uses parfor for row parallelization (requires Parallel Computing Toolbox).

seg_len = windowsize;
overlap = seg_len/2;
step = seg_len - overlap;
num_freq = floor(seg_len / 2) + 1;
frequencies = linspace(0, fs/2, num_freq);

num_rows = size(x,1);
psd = nan(num_rows, num_freq);
psd_norm = nan(num_rows, num_freq);
psd_db = nan(num_rows, num_freq);

parfor n = 1:num_rows
    x_row = x(n,:);
    validIdx = ~isnan(x_row);
    if ~any(validIdx)
        continue;
    end
    t = (0:length(x_row)-1)/fs;
    tValid = t(validIdx);
    xValid = x_row(validIdx);
    N_valid = length(xValid);
    num_segs = floor((N_valid - seg_len)/step) + 1;
    local_seg_len = seg_len;
    if num_segs < 1
        local_seg_len = N_valid;
        num_segs = 1;
    end
    sum_Pxx = zeros(num_freq, 1);
    num_used = 0;
    for s = 1:num_segs
        start_idx = 1 + (s-1)*step;
        end_idx = min(start_idx + local_seg_len - 1, N_valid);
        if end_idx - start_idx + 1 < 2, continue; end
        t_seg = tValid(start_idx:end_idx);
        x_seg = xValid(start_idx:end_idx);
        p = polyfit(t_seg, x_seg, 1);
        x_det = x_seg - polyval(p, t_seg);
        Pxx_seg = plomb(x_det, t_seg, frequencies, 'psd');
        sum_Pxx = sum_Pxx + Pxx_seg;
        num_used = num_used + 1;
    end
    if num_used == 0, continue; end
    Pxx = sum_Pxx / num_used;
    psd(n,:) = Pxx;
    var_x = var(xValid);
    if sum(Pxx) > eps
        temp_norm = Pxx / var_x;
    else
        temp_norm = Pxx;
    end
    temp_norm(temp_norm <= 0) = eps;
    psd_norm(n,:) = temp_norm;
    psd_db(n,:) = 10*log10(temp_norm);
end
if all(isnan(psd(:)))
    psd = []; frequencies = []; psd_norm = []; psd_db = [];
end
end