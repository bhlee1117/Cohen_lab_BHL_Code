function [delay yint] = get_delay(kymoTrace, binCounts, peakRange)
%[delay yint] = get_delay(kymoTrace, binCounts, peakRange)
% get_delay - Detect spike delay by spline interpolation and peak finding
% 2024.10.17 Byung Hun Lee, Cohen Lab
%
% Input:
%   kymoTrace: [N x T x S] array (ROI x Time x Trials)
%   binCounts: scalar, number of interpolation bins per sample
%   peakRange (optional): 1x2 vector [minIdx, maxIdx] for restricting peak search in original time indices
%
% Output:
%   delay: [N x S] matrix of peak timing (interpolated)

    [N, T, S] = size(kymoTrace);
    x0 = (0:T-1); % original time axis
    xint = linspace(x0(1), x0(end), binCounts*T); % interpolated axis

    if nargin < 3
        peakRange = [1, T];
    end
    % Convert to interpolated index range
    interpMin = floor((peakRange(1)-1) * binCounts) + 1;
    interpMax = ceil((peakRange(2)-1) * binCounts) + 1;

    delay = NaN(N, S); % preallocate result
yint=[];
    for s = 1:S
        for n = 1:N
            trace = kymoTrace(n,:,s);
            if ~any(isnan(trace))
                yint(n,:,s) = interp1(x0, trace, xint, 'spline');
                [~, idx] = max(yint(n,interpMin:interpMax,s));
                delay(n,s) = xint(interpMin + idx - 1);
            end
        end
    end
end
