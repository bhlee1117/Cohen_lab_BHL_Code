function [crossCorr, lags, nPairs] = nanXCorr(x, y, maxLag, normalized, minPairs)
% nanXCorr Calculates cross-correlation between two vectors while ignoring NaNs.
%
% Convention:
%   positive lag k means y is delayed relative to x:
%
%       crossCorr(k) = corr( x(t), y(t+k) )
%
% INPUT:
%   x          - First input vector, e.g. neuron i
%   y          - Second input vector, e.g. neuron j
%   maxLag     - Maximum lag in samples, optional
%   normalized - true: Pearson correlation at each lag
%                false: mean product at each lag
%   minPairs   - minimum number of valid non-NaN pairs required, optional
%
% OUTPUT:
%   crossCorr  - Cross-correlation values
%   lags       - Lags in samples
%   nPairs     - Number of valid pairs used at each lag

% Ensure column vectors
x = x(:);
y = y(:);

% Validate lengths
n = numel(x);

if numel(y) ~= n
    error('Input vectors must be of the same length.');
end

if n == 0
    error('Input vectors must not be empty.');
end

% Defaults
if nargin < 3 || isempty(maxLag)
    maxLag = n - 1;
end

if nargin < 4 || isempty(normalized)
    normalized = true;
end

if nargin < 5 || isempty(minPairs)
    minPairs = 3;
end

% Validate parameters
if ~isscalar(maxLag) || maxLag < 0 || maxLag ~= floor(maxLag)
    error('maxLag must be a nonnegative integer.');
end

if ~isscalar(minPairs) || minPairs < 1 || minPairs ~= floor(minPairs)
    error('minPairs must be a positive integer.');
end

% Avoid impossible lags
maxLag = min(maxLag, n - 1);

% Initialize outputs
lags = (-maxLag:maxLag).';
crossCorr = NaN(size(lags));
nPairs = zeros(size(lags));

for ii = 1:numel(lags)

    lag = lags(ii);

    if lag > 0
        % positive lag: y occurs after x
        % compare x(t) with y(t + lag)
        xLag = x(1:end-lag);
        yLag = y(1+lag:end);

    elseif lag < 0
        % negative lag: y occurs before x
        k = -lag;
        xLag = x(1+k:end);
        yLag = y(1:end-k);

    else
        % zero lag
        xLag = x;
        yLag = y;
    end

    % Keep only pairs where both values are observed
    validIdx = ~isnan(xLag) & ~isnan(yLag);

    xValid = xLag(validIdx);
    yValid = yLag(validIdx);

    nPairs(ii) = numel(xValid);

    % Require enough valid samples
    if nPairs(ii) < minPairs
        crossCorr(ii) = NaN;
        continue;
    end

    if normalized
        % Pearson correlation at this lag.
        % Mean subtraction is done separately for each lag/window.
        xCentered = xValid - mean(xValid);
        yCentered = yValid - mean(yValid);

        denom = sqrt(sum(xCentered.^2) * sum(yCentered.^2));

        if denom > 0
            crossCorr(ii) = sum(xCentered .* yCentered) / denom;
        else
            crossCorr(ii) = NaN;
        end

    else
        % Raw mean product at this lag, ignoring NaNs
        crossCorr(ii) = sum(xValid .* yValid) / nPairs(ii);
    end
end

end