function [psd, frequencies, psd_norm] = nanPSD(x, fs, windowSize)
    validIdx = ~isnan(x);
    xValid = x(validIdx);
    xValid = detrend(xValid);  % Remove linear trends

    if nargin < 3 || isempty(windowSize)
        windowSize = min(length(xValid), 256);  % reasonable default
    end

    win = hamming(windowSize);
    noverlap = floor(0.5 * windowSize);

    [psd, frequencies] = pwelch(xValid, win, noverlap, [], fs);
    psd_norm = psd / mean(psd) * std(x, 'omitnan')^2;
end
