function [psd, frequencies, psd_norm, psd_db] = nanPSD(x, fs, windowsize)
% nanPSD Compute PSD while handling NaNs using averaged Lomb-Scargle.
%
% [psd, frequencies, psd_norm, psd_db] = nanPSD(x, fs, windowsize)
%
% INPUT:
%   x          - Data matrix. Each row is one signal/neuron, columns are time.
%                Missing samples should be NaN.
%   fs         - Sampling frequency in Hz.
%   windowsize - Window size in samples/frames. Optional default: 3000.
%
% OUTPUT:
%   psd        - Raw PSD for each row.
%   frequencies - Frequency vector in Hz.
%   psd_norm   - Variance-normalized PSD.
%   psd_db     - 10*log10(psd_norm).
%
% Notes:
%   Windows are defined on the original regular time grid.
%   NaNs are removed only inside each window before calling plomb.
%   Progress is printed to the console every ~1% of rows completed.

% If input is a vector, treat it as one row signal
if isvector(x)
    x = x(:).';
end

[num_rows, num_time] = size(x);

% Validate fs
if nargin < 2 || isempty(fs) || ~isscalar(fs) || fs <= 0
    error('fs must be a positive scalar.');
end

% Default window size
if nargin < 3 || isempty(windowsize)
    windowsize = min(3000, num_time);
end

% Validate window size
if ~isscalar(windowsize) || windowsize < 2 || windowsize ~= floor(windowsize)
    error('windowsize must be an integer >= 2.');
end

% Do not allow window longer than the data
seg_len = min(windowsize, num_time);

overlap = floor(seg_len / 2);
step = seg_len - overlap;

% Frequency vector, similar to one-sided FFT/Welch frequencies
frequencies = (0:floor(seg_len/2)) * (fs / seg_len);
num_freq = numel(frequencies);
freq_col = frequencies(:);

% Initialize outputs
psd      = nan(num_rows, num_freq);
psd_norm = nan(num_rows, num_freq);
psd_db   = nan(num_rows, num_freq);

% Minimum number of valid points required in a window.
minValidPerWindow = 4;

% Window start positions on the original time axis
if num_time <= seg_len
    start_indices = 1;
else
    start_indices = 1:step:(num_time - seg_len + 1);
end

% Original time vector
t = (0:num_time-1) / fs;

%% ---- Progress tracking via DataQueue (compatible with parfor) ----
% parallel.pool.DataQueue lets worker threads safely send messages to the
% host without blocking. afterEach fires the callback on every send(q, ...).
q          = parallel.pool.DataQueue;
n_done     = 0;
print_step = max(1, floor(num_rows / 100));   % update every ~1%

afterEach(q, @(~) updateProgress());

    function updateProgress()
        n_done = n_done + 1;
        if mod(n_done, print_step) == 0 || n_done == num_rows
            fprintf('nanPSD: %d / %d rows complete (%.0f%%)\n', ...
                n_done, num_rows, 100 * n_done / num_rows);
        end
    end

fprintf('nanPSD: starting %d rows...\n', num_rows);

%% ---- Main loop ----
parfor n = 1:num_rows

    x_row = x(n, :);

    valid_all = isfinite(x_row);

    % Skip all-invalid rows
    if ~any(valid_all)
        send(q, n);
        continue;
    end

    sum_Pxx = zeros(num_freq, 1);
    num_used = 0;

    for s = 1:numel(start_indices)

        start_idx = start_indices(s);
        end_idx = start_idx + seg_len - 1;

        idx = start_idx:end_idx;

        x_seg_full = x_row(idx);
        valid_seg = isfinite(x_seg_full);

        % Need enough valid samples inside this time window
        if sum(valid_seg) < minValidPerWindow
            continue;
        end

        t_seg = t(idx);
        t_seg = t_seg(valid_seg);
        x_seg = x_seg_full(valid_seg);

        % Convert to column vectors
        t_seg = t_seg(:);
        x_seg = x_seg(:);

        % Use relative time for better numerical conditioning
        t_rel = t_seg - t_seg(1);

        % Linear detrend within this window
        p = polyfit(t_rel, x_seg, 1);
        x_det = x_seg - polyval(p, t_rel);

        % Lomb-Scargle PSD
        Pxx_seg = plomb(x_det, t_rel, freq_col, 'psd');
        Pxx_seg = Pxx_seg(:);

        % Skip problematic windows
        if numel(Pxx_seg) ~= num_freq || any(~isfinite(Pxx_seg))
            continue;
        end

        sum_Pxx = sum_Pxx + Pxx_seg;
        num_used = num_used + 1;
    end

    % Skip rows with no usable windows
    if num_used == 0
        send(q, n);
        continue;
    end

    % Average PSD over windows
    Pxx = sum_Pxx / num_used;

    psd(n, :) = Pxx(:).';

    % Variance-normalized PSD.
    x_valid_all = x_row(valid_all);
    var_x = var(x_valid_all);

    if isfinite(var_x) && var_x > eps
        temp_norm = Pxx / var_x;
    else
        temp_norm = nan(size(Pxx));
    end

    psd_norm(n, :) = temp_norm(:).';

    % dB version
    temp_db = temp_norm;
    temp_db(~isfinite(temp_db) | temp_db <= 0) = eps;

    psd_db(n, :) = 10 * log10(temp_db(:).');

    send(q, n);   % notify progress tracker
end

fprintf('nanPSD: done.\n');
end
