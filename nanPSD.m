function [frequencies, psd, psd_norm] = nanPSD(x, fs, windowSize)
% NANPowerSpectrum Calculates the power spectral density of a vector with NaN values.
%
% [frequencies, psd] = nanPowerSpectrum(x, fs, windowSize)
%
% INPUT:
%   x          - Input vector (can include NaN values)
%   fs         - Sampling frequency (in Hz)
%   windowSize - Window size for Welch's method (optional, default: length of x)
%
% OUTPUT:
%   frequencies - Frequency values corresponding to the PSD
%   psd         - Power spectral density values

    % Remove NaN values
    validIdx = ~isnan(x);
    xValid = x(validIdx);
    
    % Use default window size if not provided
    if nargin < 3 || isempty(windowSize)
        windowSize = length(xValid);
    end
    
    % Use Welch's method to calculate the PSD
    [psd, frequencies] = pwelch(xValid, windowSize, [], [], fs);

    %psd_norm = psd / sum(psd * (frequencies(2) - frequencies(1)));
    psd_norm= psd /(sum(psd)/length(psd)) * std(x,'omitnan')^2;
end
