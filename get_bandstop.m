function filtered_traces = get_bandstop(traces, fs, stopband)
% BANDSTOP_FILTER Applies a bandstop filter to an N x T matrix of traces.
% 
% filtered_traces = bandstop_filter(traces, fs, stopband)
% 
% Parameters:
% traces    - N x T matrix of traces, where N is the number of traces, T is the time points
% fs        - Sampling frequency (Hz)
% stopband  - 1x2 vector specifying the stopband frequencies [f_low, f_high] (Hz)
% 
% Returns:
% filtered_traces - Bandstop filtered traces (N x T)

% Validate input dimensions
if size(stopband, 2) ~= 2
    error('Stopband must be a 1x2 vector specifying [f_low, f_high].');
end

% Design the bandstop filter
f_low = stopband(1);
f_high = stopband(2);

if f_low >= f_high
    error('f_low must be less than f_high in the stopband.');
end

% Design a Butterworth bandstop filter
filter_order = 4; % Adjust as necessary for sharper or smoother filtering
[b, a] = butter(filter_order, [f_low, f_high] / (fs / 2), 'stop');

% Pre-allocate output matrix
[N, T] = size(traces);
filtered_traces = zeros(size(traces));

% Apply the filter row-wise (along each trace)
for i = 1:N
    filtered_traces(i, :) = filtfilt(b, a, traces(i, :));
end

end