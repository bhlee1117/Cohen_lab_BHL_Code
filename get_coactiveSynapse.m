function out = get_coactiveSynapse(B, n, win, varargin)
%COACTIVEAROUNDNEURON Co-activity of other neurons around neuron n's events.
%
% out = coactiveAroundNeuron(B, n, win)
% out = coactiveAroundNeuron(B, n, [pre post])
% out = coactiveAroundNeuron(..., 'ExcludeSelf', true, 'Aggregate', 'any')
%
% INPUTS
%   B   : [N x T] binary (logical or 0/1), neuron x time
%   n   : scalar neuron index (1..N)
%   win : scalar W meaning +/-W bins, OR [pre post] bins
%
% Name-Value
%   'ExcludeSelf' : default true (exclude neuron n from counts)
%   'Aggregate'   : 'any' (default) or 'sum'
%       'any' = count neuron if it fired at least once in window
%       'sum' = total number of spikes across all other neurons in window
%
% OUTPUTS (struct)
%   out.eventTimes         : times where neuron n is active
%   out.countPerEvent      : [K x 1] count per event time
%   out.fractionPerEvent   : [K x 1] count/(N-1 or N)
%   out.countPerTime       : [1 x T] count for every time bin (NaN where n inactive)
%   out.window             : [pre post]
%
% Notes: If neuron n has consecutive 1s, each bin is treated as an event.
%        If you want event onsets only, set Bn = diff([0 B(n,:)])==1 before.

p = inputParser;
p.addParameter('ExcludeSelf', true, @(x)islogical(x) && isscalar(x));
p.addParameter('Aggregate', 'any', @(s)ischar(s) || isstring(s));
p.parse(varargin{:});
excludeSelf = p.Results.ExcludeSelf;
agg = lower(string(p.Results.Aggregate));

B = logical(B);
[N,T] = size(B);

if isscalar(win)
    pre = win; post = win;
else
    pre = win(1); post = win(2);
end
pre = max(0, round(pre));
post = max(0, round(post));

Bn = B(n,:);
eventTimes = find(Bn);

K = numel(eventTimes);
countPerEvent = zeros(K,1);

% Choose which neurons to consider
idxNeurons = true(N,1);
if excludeSelf, idxNeurons(n) = false; end
M = sum(idxNeurons);

% Pre-extract other neurons
Other = B(idxNeurons, :);

for i = 1:K
    t = eventTimes(i);
    a = max(1, t-pre);
    b = min(T, t+post);

    W = Other(:, a:b);

    switch agg
        case "any"
            countPerEvent(i) = sum(any(W, 2));     % neurons co-active at least once
        case "sum"
            countPerEvent(i) = sum(W(:));          % total spikes across others
        otherwise
            error("Aggregate must be 'any' or 'sum'.");
    end
end

fractionPerEvent = countPerEvent / max(M,1);

countPerTime = nan(1,T);
countPerTime(eventTimes) = countPerEvent;

out = struct();
out.eventTimes = eventTimes;
out.countPerEvent = countPerEvent;
out.fractionPerEvent = fractionPerEvent;
out.countPerTime = countPerTime;
out.window = [pre post];
out.n = n;
out.aggregate = char(agg);
end