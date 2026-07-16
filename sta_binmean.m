function [bs, STA] = sta_binmean(dFF, trig, woff, OH, cnt, shift)
% STA_BINMEAN  Trigger-averaged trace, binned in time, memory-capped.
%   [bs, STA] = sta_binmean(dFF, trig, woff, OH, cnt, shift)
%   dFF   : nROI x nT   traces
%   trig  : 1 x K       trigger frames (must satisfy trig+woff within [1,nT] for shift=0)
%   woff  : 1 x winLen  frame offsets around each trigger (e.g. -pre:post)
%   OH    : winLen x nBin sparse one-hot mapping STA frames -> time bins
%   cnt   : 1 x nBin    # STA frames per bin (0-bins set to NaN by caller)
%   shift : scalar circular shift of the TRACE (applied as a wrapped trigger offset)
%
%   STA : nROI x winLen  trigger-average trace
%   bs  : nROI x nBin    STA averaged within each time bin
% Windows are gathered in trigger-batches so peak memory is nROI*winLen*batch,
% not nROI*K*winLen (which get_STA would allocate).
[Nsyn, nT] = size(dFF);
winLen = numel(woff);  woff = woff(:)';
K = numel(trig);
if K==0, STA = nan(Nsyn,winLen); bs = nan(Nsyn,size(OH,2)); return; end
sumMat = zeros(Nsyn, winLen);
bsz = 256;
for a = 1:bsz:K
    tb  = trig(a:min(K,a+bsz-1)); tb = tb(:);      % nb x 1
    idx = mod((tb + woff) - 1 - shift, nT) + 1;    % nb x winLen, wrapped
    chunk = reshape(dFF(:, idx'), Nsyn, winLen, []);  % Nsyn x winLen x nb
    sumMat = sumMat + sum(chunk, 3);
end
STA = sumMat / K;
bs  = (STA * OH) ./ cnt;
end
