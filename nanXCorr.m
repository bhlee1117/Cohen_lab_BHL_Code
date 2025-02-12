function [crossCorr, lags] = nanXCorr(x, y, maxLag, normalized)
% NANCrossCorr Calculates the cross-correlation between two vectors, optionally normalized, ignoring NaN values.
%
% [crossCorr, lags] = nanCrossCorr(x, y, maxLag, normalized)
%
% INPUT:
%   x          - First input vector
%   y          - Second input vector
%   maxLag     - Maximum lag to consider (optional, default: length of the vectors - 1)
%   normalized - Boolean to indicate whether to normalize (optional, default: true)
%
% OUTPUT:
%   crossCorr - Cross-correlation values
%   lags      - Lags corresponding to the cross-correlation values

% Ensure inputs are column vectors
x = x(:);
y = y(:);

% Validate lengths
if length(x) ~= length(y)
    error('Input vectors must be of the same length.');
end

% Default maxLag if not provided
if nargin < 3
    maxLag = length(x) - 1;
end

% Default normalization if not provided
if nargin < 4
    normalized = true;
end

% Initialize output
lags = -maxLag:maxLag;
crossCorr = zeros(size(lags));

% Calculate cross-correlation
for i = 1:length(lags)
    lag = lags(i);

    if lag < 0
        % Shift x forward
        validIdx = ~isnan(x(1:end+lag)) & ~isnan(y(-lag+1:end));
        xValid = x(1:end+lag);
        yValid = y(-lag+1:end);
        if isempty(find(validIdx, 1))
            crossCorr(i) = NaN;
        else
            if normalized
                crossCorr(i) = corr(xValid(validIdx), yValid(validIdx));
            else
                crossCorr(i) = sum(xValid(validIdx) .* yValid(validIdx))/length(validIdx);
            end
        end
    elseif lag > 0
        % Shift y forward
        validIdx = ~isnan(x(lag+1:end)) & ~isnan(y(1:end-lag));
        xValid = x(lag+1:end);
        yValid = y(1:end-lag);
        if isempty(find(validIdx, 1))
            crossCorr(i) = NaN;
        else
            if normalized
                crossCorr(i) = corr(xValid(validIdx), yValid(validIdx));
            else
                crossCorr(i) = sum(xValid(validIdx) .* yValid(validIdx))/length(validIdx);
            end
        end
    else
        % No shift
        validIdx = ~isnan(x) & ~isnan(y);
        if isempty(find(validIdx, 1))
            crossCorr(i) = NaN;
        else
            if normalized
                crossCorr(i) = corr(x(validIdx), y(validIdx));
            else
                crossCorr(i) = sum(x(validIdx) .* y(validIdx))/length(validIdx);
            end
        end
    end
end
end
