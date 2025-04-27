function filtered_traces = get_bandpass(traces, fs, passband)
% GET_BANDPASS Applies a bandpass filter to an N x T matrix of traces.
% Handles NaNs by filtering only non-NaN segments.
% 
% filtered_traces = get_bandpass(traces, fs, passband)
% 
% Parameters:
% traces    - N x T matrix of traces (double), may contain NaNs
% fs        - Sampling frequency (Hz)
% passband  - 1x2 vector specifying the passband frequencies [f_low, f_high] (Hz)
% 
% Returns:
% filtered_traces - Bandpass filtered traces (N x T), NaNs are preserved

% Validate passband input
if size(passband, 2) ~= 2
    error('Passband must be a 1x2 vector specifying [f_low, f_high].');
end

f_low = passband(1);
f_high = passband(2);

if f_low >= f_high
    error('f_low must be less than f_high in the passband.');
end

% Design the Butterworth bandpass filter
filter_order = 4;
[b, a] = butter(filter_order, [f_low, f_high] / (fs / 2), 'bandpass');

[N, T] = size(traces);
filtered_traces = NaN(N, T);  % Initialize with NaNs to preserve them

for i = 1:N
    trace = traces(i, :);
    
    % Find valid (non-NaN) segments
    isnan_mask = isnan(trace);
    valid_mask = ~isnan_mask;
    
    % Identify contiguous non-NaN segments
    d = diff([false, valid_mask, false]);
    starts = find(d == 1);
    ends = find(d == -1) - 1;

    for j = 1:length(starts)
        seg_start = starts(j);
        seg_end = ends(j);
        segment = trace(seg_start:seg_end);
        
        if length(segment) > 6 * filter_order  % Require enough length for filtfilt
            filtered_segment = filtfilt(b, a, segment);
            filtered_traces(i, seg_start:seg_end) = filtered_segment;
        else
            % If too short to filter, leave as NaN or copy original values
            filtered_traces(i, seg_start:seg_end) = NaN;
        end
    end
end

end
