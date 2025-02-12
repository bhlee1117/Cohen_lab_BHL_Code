function [wt f_int]= get_waveletTransform(signal,fs,freqRange)

if nargin < 3
    freqRange = [0 fs/2];
    f_int = 0:0.1:fs/2;
else
    f_int = freqRange(1):0.1:freqRange(2);
end

nan_indices = isnan(signal);
segments = bwlabel(~nan_indices); % Group contiguous non-NaN segments

wt = []; % Initialize wavelet transform matrix

% Threshold for short segments
short_segment_threshold = 40; % Minimum length for CWT

% Process each segment
for i = 1:max(segments)
    segment = signal(segments == i);
    
    if ~isempty(segment)
        if length(segment) >= short_segment_threshold
            % Use CWT for long segments
            [wt_tmp, f] = cwt(segment, 'amor', fs, 'FrequencyLimits', freqRange);
            wt(:, find(segments == i)) = interp1(f, wt_tmp, f_int, 'linear', 'extrap');
        else
            % Mark frames corresponding to short segments as NaN
            wt(:, find(segments == i)) = NaN;
        end
    end
end

wt(:, nan_indices) = NaN; % Reinsert NaN values where the input signal had NaNs

end
