function filtered_traces = get_bandstop(traces, fs, stopband)
%filtered_traces = get_bandstop(traces, fs, stopband)
% GET_BANDSTOP Applies a bandstop filter to an N x T matrix of traces.
% Handles NaNs by filtering only non-NaN segments.
%
% filtered_traces = get_bandstop(traces, fs, stopband)
%
% Parameters:
% traces    - N x T matrix of traces (double), may contain NaNs
% fs        - Sampling frequency (Hz)
% stopband  - 1x2 vector specifying the stopband frequencies [f_low, f_high] (Hz)
%
% Returns:
% filtered_traces - Bandstop filtered traces (N x T), NaNs are preserved

% Validate stopband input
if size(stopband, 2) ~= 2
    error('Stopband must be a 1x2 vector specifying [f_low, f_high].');
end

f_low = stopband(1);
f_high = stopband(2);

if f_low >= f_high
    error('f_low must be less than f_high in the stopband.');
end

% Design Butterworth bandstop filter
filter_order = 4;
[b, a] = butter(filter_order, [f_low, f_high] / (fs / 2), 'stop');

[N, T] = size(traces);
filtered_traces = NaN(N, T);  % Initialize with NaNs to preserve them

for i = 1:N
    trace = traces(i, :);

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

        if length(segment) > 3 * filter_order
            filtered_segment = filtfilt(b, a, segment);
            filtered_traces(i, seg_start:seg_end) = filtered_segment;
        else
            filtered_traces(i, seg_start:seg_end) = NaN;
        end
    end
end

end
