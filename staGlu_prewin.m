function [STA, tMs, score, pval] = staGlu_prewin(dFF, trigIdx, dtG, winMs, scoreMs, nShuffle, seed)
% STAGLU_PREWIN  Spike-triggered glutamate STA + windowed score(s) with shuffle p-value.
%   [STA, tMs, score, pval] = staGlu_prewin(dFF, trigIdx, dtG, winMs, scoreMs, nShuffle, seed)
%   dFF      : Nsyn x nTg  glutamate traces
%   trigIdx  : trigger glutamate frames (e.g. V-spikes mapped to glu frames)
%   dtG      : s / glutamate frame
%   winMs    : [pre post] ms for the STA window (both positive)
%   scoreMs  : W x 2 matrix of [a b] ms windows (relative to spike) scored for tuning,
%              e.g. [-50 0] (pre-spike) or [-50 0; 0 50] (pre- AND post-spike). Each row
%              is one window; results are returned column-wise in the same order.
%   nShuffle : # circular shuffles                seed : rng seed (default 0)
%
%   STA   : Nsyn x winLen  trigger-averaged dF/F        tMs : 1 x winLen  (ms from spike)
%   score : Nsyn x W  mean dF/F in each scoreMs window (observed)
%   pval  : Nsyn x W  right-tailed circular-shuffle p-value for each column of `score`
% Circular shuffle rolls the glutamate trace (triggers fixed): tests whether the
% glutamate in each window is elevated above chance. The same null (one set of shifts)
% is scored for every window, so the shuffle is computed only once.
if nargin<7, seed = 0; end
nTg = size(dFF,2);
pre = round(winMs(1)/1000/dtG); post = round(winMs(2)/1000/dtG);
woff = -pre:post; winLen = numel(woff); tMs = woff*dtG*1000;

W = size(scoreMs,1);                                              % # of score windows
rIdx = []; cIdx = []; cnt = zeros(1,W);
for w = 1:W
    inW    = find(tMs>=scoreMs(w,1) & tMs<=scoreMs(w,2));
    rIdx   = [rIdx; inW(:)];   cIdx = [cIdx; w*ones(numel(inW),1)]; %#ok<AGROW>
    cnt(w) = numel(inW);
end
OH = sparse(rIdx, cIdx, 1, winLen, W);       % winLen x W one-hot -> per-window mean
cnt(cnt==0) = NaN;

trigIdx = trigIdx(trigIdx+woff(1)>=1 & trigIdx+woff(end)<=nTg);    % keep in-range windows
if isempty(trigIdx)
    STA = nan(size(dFF,1),winLen); score = nan(size(dFF,1),W); pval = nan(size(dFF,1),W); return;
end

[score, STA] = sta_binmean(dFF, trigIdx, woff, OH, cnt, 0);       % score: Nsyn x W

minShift = max(1, round(1/dtG)); minShift = min(minShift, max(1,floor(nTg/2)));
hiShift  = max(minShift+1, nTg-minShift);
rng(seed); shifts = randi([minShift, hiShift], 1, nShuffle);
nullS = zeros(size(dFF,1), W, nShuffle);
for i = 1:nShuffle
    nullS(:,:,i) = sta_binmean(dFF, trigIdx, woff, OH, cnt, shifts(i));
end
pval = (1 + sum(nullS >= score, 3)) ./ (nShuffle + 1);            % Nsyn x W
pval(isnan(score)) = NaN;
end
