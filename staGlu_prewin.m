function [STA, tMs, score, pval] = staGlu_prewin(dFF, trigIdx, dtG, winMs, scoreMs, nShuffle, seed)
% STAGLU_PREWIN  Spike-triggered glutamate STA + pre-window score with shuffle p-value.
%   [STA, tMs, score, pval] = staGlu_prewin(dFF, trigIdx, dtG, winMs, scoreMs, nShuffle, seed)
%   dFF      : Nsyn x nTg  glutamate traces
%   trigIdx  : trigger glutamate frames (e.g. V-spikes mapped to glu frames)
%   dtG      : s / glutamate frame
%   winMs    : [pre post] ms for the STA window (both positive)
%   scoreMs  : [a b] ms window (relative to spike) used as the tuning score, e.g. [-50 0]
%   nShuffle : # circular shuffles                seed : rng seed (default 0)
%
%   STA   : Nsyn x winLen  trigger-averaged dF/F        tMs : 1 x winLen  (ms from spike)
%   score : Nsyn x 1  mean dF/F in scoreMs window (observed)
%   pval  : Nsyn x 1  right-tailed circular-shuffle p-value for `score`
% Circular shuffle rolls the glutamate trace (triggers fixed): tests whether the
% pre-spike glutamate is elevated above chance.
if nargin<7, seed = 0; end
nTg = size(dFF,2);
pre = round(winMs(1)/1000/dtG); post = round(winMs(2)/1000/dtG);
woff = -pre:post; winLen = numel(woff); tMs = woff*dtG*1000;

inScore = tMs>=scoreMs(1) & tMs<=scoreMs(2);
OH  = sparse(find(inScore), ones(nnz(inScore),1), 1, winLen, 1);   % winLen x 1 -> pre-window mean
cnt = nnz(inScore); if cnt==0, cnt = NaN; end

trigIdx = trigIdx(trigIdx+woff(1)>=1 & trigIdx+woff(end)<=nTg);    % keep in-range windows
if isempty(trigIdx)
    STA = nan(size(dFF,1),winLen); score = nan(size(dFF,1),1); pval = nan(size(dFF,1),1); return;
end

[score, STA] = sta_binmean(dFF, trigIdx, woff, OH, cnt, 0);

minShift = max(1, round(1/dtG)); minShift = min(minShift, max(1,floor(nTg/2)));
hiShift  = max(minShift+1, nTg-minShift);
rng(seed); shifts = randi([minShift, hiShift], 1, nShuffle);
nullS = zeros(size(dFF,1), nShuffle);
for i = 1:nShuffle
    nullS(:,i) = sta_binmean(dFF, trigIdx, woff, OH, cnt, shifts(i));
end
pval = (1 + sum(nullS >= score, 2)) ./ (nShuffle + 1);
pval(isnan(score)) = NaN;
end
